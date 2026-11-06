#!/usr/bin/env python3
"""Linux regression for the vmatch inherited-descriptor launcher."""
import os
import re
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

if not sys.platform.startswith("linux") or not Path("/proc/self/fd").is_dir():
    print("vmatch descriptor regression skipped: Linux procfs probes required")
    raise SystemExit(0)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = (ROOT / "CRISPRCasFinder.pl").read_text()


def subroutine(text, name):
    start = text.index("sub " + name + "\n")
    end = text.index("\nsub ", start + 5)
    return text[start:end]


new = SOURCE[SOURCE.index("sub makevmatchcall\n"):SOURCE.index("\nsub makesystemcall\n")]
old = subroutine(SOURCE, "makesystemcall")
assert SOURCE.count('makevmatchcall("vmatch "') == 1
assert 'makesystemcall("vmatch "' not in SOURCE

DRIVER = """use strict; use warnings;
{old}
{new}
my ($mode, $command, $check_parent) = @ARGV;
if ($mode eq 'new') {{ makevmatchcall($command); }}
else {{ makesystemcall($command); }}
die "parent fd missing\\n" if $check_parent && !-e "/proc/$$/fd/68";
print "parent-ok\\n";
""".format(old=old, new=new)


def run(perl, script, mode, command, fds=(), cwd=None):
    args = [perl, str(script), mode, command] + (["check"] if fds else [])
    return subprocess.run(args, pass_fds=tuple(fds),
                          capture_output=True, text=True, timeout=10, cwd=cwd)


with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    driver = td / "driver.pl"
    driver.write_text(DRIVER)
    perl = os.environ.get("PERL", "perl")
    # The driver requires Perl; PERL permits maintainers to select one.
    if subprocess.run([perl, "-e", "1"], capture_output=True, timeout=10).returncode:
        raise SystemExit("selected PERL cannot run")

    handles = [os.open(os.devnull, os.O_RDWR) for _ in range(66)]
    try:
        assert 68 in handles, "run standalone so descriptor 68 is among the test handles"
        result = run(perl, driver, "new", "test ! -e /proc/$$/fd/68", handles)
        assert result.returncode == 0, result.stderr
        assert result.stdout == "parent-ok\n"
        result = run(perl, driver, "new", "sh -c 'exit 7'", handles)
        assert result.returncode == 1 and "errorcode 1792" in result.stderr
        result = run(perl, driver, "new", "kill -TERM $$", handles)
        assert result.returncode == 1 and re.search(r"errorcode 15(?:\n|$)", result.stderr)
    finally:
        for fd in handles:
            os.close(fd)
    print("core descriptor, parent-handle and failure checks passed")

    native = os.environ.get("VMATCH")
    mkvtree = os.environ.get("MKVTREE")
    if not native or not mkvtree:
        print("native vmatch regression skipped: set VMATCH and MKVTREE")
        raise SystemExit(0)
    fixture = td / "fixture.fna"
    fixture.write_text(">synthetic\n" + "ACGT" * 2000 + "\n")
    subprocess.run([mkvtree, "-db", str(fixture), "-dna", "-pl", "-lcp", "-suf",
                    "-tis", "-ois", "-bwt", "-bck", "-sti1"], check=True,
                   cwd=td, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10)
    command = "{} -l 3 -s leftseq -evalue 1 -absolute -nodist -noevalue -noscore -noidentity -sort ia -best 10 {} > result.txt".format(shlex.quote(native), shlex.quote(str(fixture)))
    clean = run(perl, driver, "new", command, cwd=td)
    assert clean.returncode == 0, clean.stderr
    output = (td / "result.txt").read_bytes()
    assert output
    handles = [os.open(os.devnull, os.O_RDWR) for _ in range(66)]
    try:
        assert 68 in handles, "run standalone so descriptor 68 is among the test handles"
        baseline = run(perl, driver, "old", command, handles, cwd=td)
        failure = re.search(r"filedescriptor \d+ is too large", baseline.stderr)
        assert baseline.returncode != 0 and failure, baseline.stderr
        print("unmodified launch reproduced: " + failure.group(0))
        fixed = run(perl, driver, "new", command, handles, cwd=td)
        assert fixed.returncode == 0 and (td / "result.txt").read_bytes() == output
    finally:
        for fd in handles:
            os.close(fd)
    print("native vmatch regression passed")

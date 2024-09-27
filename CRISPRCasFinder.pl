#!/usr/bin/env perl

######################################################################################
#  CRISPRCasFinder, an update of CRISRFinder, includes a portable version,
#  enhanced performance and integrates search for cas genes.
# 
#  Copyright (C) 2017- CRISPR-Cas++ team 
#  (CNRS, Université Paris-Saclay I2BC, Institut Pasteur C3BI, Université de Lille)
#
#  See the COPYRIGHT file for details.
#
#  CRISPRCasFinder is distributed under the terms of the GNU General Public License
#  (GPLv3). See the COPYING file for details.  
######################################################################################

# CPAN modules
## insert perlib ##
use strict;
use Class::Struct;   #charge le module qui construit les struct
use warnings;
use File::Copy;
use File::Basename;
use Data::Dumper;
use Cwd;
#use Date::Calc qw(:all);
use Getopt::Long;
#use Unix::Sysexits;

# Bio Perl modules
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::DB::Fasta; #to extract sequence from fasta file

# TODO add URLescape module in case the GFF attributes get messy
#local modules
my $version = "4.3.2";

# set parameters for the Vmatch program
#$ENV{'LD_LIBRARY_PATH'} = '.';

## Parameters
my $SpSim = 60; # maximal allowed percentage of similarity between Spacers (default value=60)

my $Sp2 = 2.5; # maximal Spacers size in function of DR size (default value=2.5)

my $Sp1 = 0.6; # minimal Spacers size in function of DR size (default value=0.6)

my $mismOne = 1; # allow mismatchs (default value=1)

my $S2 = 60; # maximal size of Spacers (default value=60)

my $S1 = 25; # minimal size of Spacers (default value=25)

my $M2 = 55; # maximal size of DRs (default value=55)

my $M1 = 23; # minimal size of DRs (default value=23)

my $DRtrunMism = 33.3; # %mismatchs allowed for truncated DR (default value=33.3)

my $DRerrors = 20; # %mismatchs allowed between DRs (default value=20)

my $userfile = "";

my $launchCasFinder = 0; # boolean variable indicating if we use casfinder or not (default value=0)

#my $casfinder = "CasFinder-2.0.3"; # repository containing new CasFinder (default 'CasFinder-2.0.3')

my $kingdom = "Bacteria"; # allow to choose analysis between Archaea and Bacteria (default 'Bacteria')

my $vicinity = 600; # number of nucleotides separating CRISPR array from next Cas system (default value=600)

my $writeFullReport = 0; # boolean variable indicating if we write crispr-cas_vicinity_report or not (default value=0)

my $so = "./sel392v2.so"; # path to shared object (.so) file (former name: $pathSoFile)

my $crisprdb = ""; # path to all CRISPR candidates contained in CRISPRdb (from last update)

my $repeats = ""; # path to file containing repeat sequences, IDs, and Nb in CRISPRdb (last update)

my $dirRepeat = ""; # path to file containing repeat IDs and Orientation according to CRISPRDirection

my $html = 0; # boolean variable indicating if we use html visualization or not (default value=0)

my $outputDirName = ""; # repository name containing results to be set by user (default: word 'Result_' followed by basename of input file)

my $flankingRegion = 100; # size of flanking regions in base pairs (bp) for each analyzed CRISPR array (default value=100)

my $cpuProkka = 1; # number of CPUs to use for Prokka (default value=1)

my $cpuMacSyFinder = 1; # number of CPUs to use for MacSyFinder (default value=1)

my $rcfowce = 0; # option allowing to run CasFinder only when any CRISPR exists (default value=0) (set if -cas is set)

my $metagenome = 0; # option allowing to analyze metagenome (default value=0)

my $logOption = 0; # option allowing to write LOG files (default value=0)

my $keep = 0; # option allowing to keep secondary folders/files (Prodigal/Prokka, CasFinder, rawFASTA, Properties); default value=0

my $definition = "SubTyping"; # option allowing to specify CasFinder definition (if option -cas is set) to be more or less stringent (default value="SubTyping"). Other options include "Typing" and "General". 

my $userGFF = ""; # option allowing user to provide an annotation GFF file (if options -cas and -faa are set) (default value='')

my $userFAA = ""; # option allowing user to provide a proteome file '.faa' (if options -cas and -gff are set) (default value='')

my $clusteringThreshold = 0; # option allowing (if option -cas is set) to constitute clusters or groups of CRISPR or Cas systemes given a determined threshold e.g. 20kb (default value=0)

my $useProdigal = 1; # option allowing to use Prodigal as default option instead of Prokka (default value=1)

my $useProkka = 0; # option allowing to use Prokka instead of Prodigal (default value=0)

my $gscf = 0; # option allowing to get summary file of Cas-finder and copy it to TSV repository (default value=0) 

my $cssFile = ""; # option allowing to copy CSS file (crispr.css) to get the same design as CRISPRdb when using option -HTML (default value='')

my $genCode = 11; # option allowing to modify the genetic code (translation table) for CDS annotation (default value=11)

my $levelMin = 1; # option allowing to choose the minimum evidence-level corresponding to CRISPR arrays we want to display (default value=1)

my $quiet = 0; # option allowing to run the program quieter (default value=0) 

my $seqMinSize = 0; #NV option allowing to fix a sequence minimum size to search CRISPR and Cas systems (lighter process on big Data) (default value=0)

my $fast = 0; # option allowing to run the program faster (default value=0)

my $force = 0; # option allowing to force/foster detection of specific CRISPR arrays (default value=0)

my $fosteredDRLength = 30; #option allowing to foster a specific repeat's length when option '-force' is set (default value=30)

my $fosteredDRBegin = "G"; #option allowing to foster a specific repeat's beginning when option '-force' is set (default value='G'), regular expressions could be considered

my $fosteredDREnd = "AA."; #option allowing to foster a specific repeat's ending when option '-force' is set (default value='AA.'), regular expressions could be considered 

my $repeatsQuery = ""; # option allowing to specify a query file containing repeats to be matched (default value='')

my $minNbSpacers = 1; # option allowing to specifiy the minimum number of spacers required for each array (default value=1)

my $betterDetectTruncatedDR = 0; # option allowing to better detect the truncated DR (default value=0)

my $percentageMismatchesHalfDR = 4; # option allowing to set the percentage of allowed mismatches in truncated DR (default value=4)

my $classifySmall = 0; # option allowing to change evidence level status of small arrays (evidence-level 1) having the same consensus repeat as an evidence-level 4 array (default value=0)

# Variables used to determine mode
my $print_version = 0;  # Print version
my $print_help = 0;     # Print help

##

## Manage arguments
# Note that standard overloading is possible i.e.
# using -cpuMacSysFinder after -fast will set the
# option to the desired value (it was the case with the
# previous implementation).
# Note that $DRtrunMism and $DRerrors are modified after.
GetOptions (
    # General options
    "help|h" => \$print_help,
    "version|v" => \$print_version,
    # Input/Output and -so
    "in|i=s" => \$userfile,
    "outdir|out=s" => \$outputDirName,
    "keepAll|keep" => \$keep,
    "LOG|log" => \$logOption,
    "HTML|html" => \$html,
    "copyCSS=s" => \$cssFile,
    "soFile|so=s" => \$so,
    "quiet|q" => \$quiet,
    "faster|fast" => sub {
	# Fast mode
	$fast = 1;
	# Use all CPU for Prokka and MacSyFinder
	$cpuProkka = 0;
	$cpuMacSyFinder = 0;
    },
    "minSeqSize|mSS=i" => \$seqMinSize,
    # Options for detection of CRISPR arrays
    "mismDRs|md=f" => \$DRerrors,
    "truncDR|t=f" => \$DRtrunMism,
    "minDR|mr=i" => \$M1,
    "maxDR|xr=i" => \$M2,
    "minSP|ms=i" => \$S1,
    "maxSP|xs=i" => \$S2,
    "noMism|n" => sub { $mismOne=0; },
    "percSPmin|pm=f" => \$Sp1,
    "percSPmax|px=f" => \$Sp2,
    "spSim|s=f" => \$SpSim,
    "DBcrispr|dbc=s" => \$crisprdb,
    "repeats|rpts=s" => \$repeats,
    "DIRrepeat|drpt=s" => \$dirRepeat,
    "flank|fl=i" => \$flankingRegion,
    "levelMin|lMin=i" => \$levelMin,
    "classifySmallArrays|classifySmall|cSA" => \$classifySmall,
    "forceDetection|force" => sub {
	$force = 1;
	$M1 = $fosteredDRLength;
    },
    "fosterDRLength|fDRL=i" => \$fosteredDRLength,
    "fosterDRBegin|fDRB=s" => \$fosteredDRBegin,
    "fosterDREnd|fDRE=s" => \$fosteredDREnd,
    "MatchingRepeats|Mrpts=s" => \$repeatsQuery,
    "minNbSpacers|mNS=i" => \$minNbSpacers,
    "betterDetectTrunc|bDT" => \$betterDetectTruncatedDR,
    "PercMismTrunc|PMT=f" => \$percentageMismatchesHalfDR,
    # Options for detection of Cas clusters
    "cas|cs" => \$launchCasFinder,
    "ccvRep|ccvr" => \$writeFullReport,
    "vicinity|vi=i" => \$vicinity,
    #"CASFinder|cf|CasFinder=s" => \$casfinder,
    "cpuMacSyFinder|cpuM=i" => \$cpuMacSyFinder,
    "rcfowce" => \$rcfowce,
    "definition|def=s" => \$definition, # Values: "SubTyping", "Typing", "General"
    "gffAnnot|gff=s" => \$userGFF,
    "proteome|faa=s" => \$userFAA, # Need -cas and -gff
    "cluster|ccc=i" => \$clusteringThreshold, # Need -cas
    "getSummaryCasfinder|gscf" => \$gscf,
    "geneticCode|gcode=i" => \$genCode,
    "metagenome|meta" => \$metagenome,
    # Options to use Prokka
    "useProkka|prokka" => sub { $useProkka = 1; $useProdigal = 0; },
    "cpuProkka|cpuP=i" => \$cpuProkka,
    "ArchaCas|ac" => sub { $launchCasFinder=1; $kingdom = "Archaea"; },
) or do {
    print STDERR "Error in command line arguments\n";
    printhelpbasic($0);
    exit(1);
};

# If a positional argument remains, it's the input file
my $n_rem_arguments = @ARGV;
if ($n_rem_arguments) {
    $userfile = $ARGV[0];
}

# Print help and exit
if ($print_help)
{
    printhelpall($0);
    exit(0);
}

# Print version and exit
if ($print_version)
{
    printversion($0);
    exit(0);
}

# Basic checks of CLI arguments
# Output dir is created later if not given

# Check that input file was given and exists
if (!$userfile) {
    print STDERR "Error: input file not given.\n";
    printhelpbasic($0);
    exit(1);
}
if (not -e $userfile) {
    print STDERR "Error: file $userfile not found. Please check that the file exists or enter a correct file name.\n";
    printhelpbasic($0);
    exit(1);
}

print "################################################################\n";
print "# --> Welcome to $0 (version $version) \n";
print "################################################################\n\n\n";

## Control dependencies

my $vmatchProg = isProgInstalled("vmatch");
if($vmatchProg){
  print "vmatch is...............OK \n";
}
else
{
  print "vmatch is not installed, please install it and try again.\n";
  exit(1);
}

my $mkvtreeProg = isProgInstalled("mkvtree");
if($mkvtreeProg){
  print "mkvtree is...............OK \n";
}
else
{
  print "mkvtree (Vmatch dependency) is not installed, please install it and try again.\n";
  exit(1);
}

my $vsubseqselectProg = isProgInstalled("vsubseqselect");
if($vsubseqselectProg){
  print "vsubseqselect is...............OK \n";
}
else
{
  print "vsubseqselect (Vmatch dependency) is not installed, please install it and try again.\n";
  exit(1);
}
#--
my $fuzznuc = isProgInstalled("fuzznuc");
if($fuzznuc){
  print "fuzznuc (from emboss) is...............OK \n";
}
else
{
    print "fuzznuc (from emboss) is not installed, please install it and try again.\n";
    exit(1);
}

my $needle = isProgInstalled("needle");
if($needle){
  print "needle (from emboss) is...............OK \n\n\n";
}
else
{
  print "needle (from emboss) is not installed, please install it and try again.\n";
  exit(1);
}

##


## Manage DRs
$DRtrunMism = 100/$DRtrunMism;
$DRerrors = $DRerrors/100;
##

## Time - Begin
#my ($start_hour,$start_min,$start_sec) = Now(); # Now([+1]) replaced by Now()
#my($start_year,$start_month,$start_day, $start_hour,$start_min,$start_sec) = Today_and_Now();
##
my $epoch = time();  # start time

my $seqIO = Bio::SeqIO->new(-format=>'Fasta', -file=>$userfile);
my $inputfileCount=0;
my ($seq,$inputfile);
my $basename = basename($userfile);
my @outdir = split /\./, $basename;
my $outdir = $outdir[0];
#$outdir .= "_".$start_day."_".$start_month."_".$start_year."_".$start_hour."_".$start_min."_".$start_sec;
$outdir .= "_".$epoch;

# DC - 11/05/2017 Retrieve date (localtime) and default repository
my $ResultDir = "";

if($outputDirName eq ""){
  $ResultDir = "Result_".$outdir;
}
else{
  $ResultDir = $outputDirName;
}

my $ResultDirFinal = $ResultDir;

# LK
if(-d $ResultDirFinal){
  die "Output directory $ResultDirFinal already exists. Remove the directory before continuing.\nSuggestion:\nrm -rf $ResultDirFinal/\n";
}

# $ResultDir now corresponds to the temporary result directory ($outdir)
#$ResultDir = $outdir;
mkdir $ResultDir unless -d $ResultDir;


mkdir $ResultDir."/CRISPRFinderProperties" unless -d $ResultDir."/CRISPRFinderProperties"; #NV


#create a directory GFF and move all GFFs to this directory
mkdir $ResultDir."/GFF" unless -d $ResultDir."/GFF";

struct Rep => {
    Pos1  => '$',
    Length  => '$',
    DRseq => '$',
};


my $htmlFile = "index.html"; # web page allowing a simple visualization of CRISPRs and Cas genes ($ResultDir/ removed)
if($html){
  my @status = stat($userfile);
  #open (HTML, ">$htmlFile") or die "open : $!";
  open(HTML,">",$htmlFile) or die("Could not open the file $htmlFile because $!\n");
  
  print HTML "<!DOCTYPE html>\n";
  print HTML "<html>\n";
  print HTML "<head>\n";
  print HTML "<title>CRISPR-Cas viewer</title>\n";
  print HTML "<link rel='stylesheet' type='text/css' href='crispr.css'>\n";
  print HTML "</head>\n";
  print HTML "<body>\n";

  print HTML "<div class = content>\n";
  print HTML "<h1> CRISPRs arrays and Cas genes found in the submitted sequence(s) </h1>\n";
  print HTML "<br/>\n";
  print HTML "<pre>\n";
  print HTML "File name           : <font COLOR= \#FF0000> $userfile </font>\n";
  print HTML "File size (in bytes): $status[7]\n";
  print HTML "</pre>\n";

}

my @statusLOG = stat($userfile); # Get input file size in bytes

# Two additional log files (logfile and logSequences) -DC 06/2017  
#my $logfile = "logFile_".$start_day."_".$start_month."_".$start_year."_".$start_hour."_".$start_min."_".$start_sec.".txt"; # to write command lines
#my $logSeq = "logSequences_".$start_day."_".$start_month."_".$start_year."_".$start_hour."_".$start_min."_".$start_sec.".tsv"; # to write quick information on sequences
my $logfile = "logFile_".$epoch;
my $logSeq = "logSequences_".$epoch;

if($logOption){
	open (LOG, ">$logfile") or die "open : $!";
	open (LOGSEQ, ">$logSeq") or die "open : $!";
	print LOGSEQ "Sequence ID\tName\tSize (bp)\tAT%\tNb_CRISPRs\tTotal Time (s)\tCRISPR Time (s)\tCas Time (s)\tFile size (bytes)\n";
}

#my $jsonSeq = "sequences.json"; # JSON file to get information on sequence name and sequence length (bp)
#open (JSONSEQ, ">$jsonSeq") or die "open : $!";
#print JSONSEQ "[\n";
#my $jsonLineSeq = "";

## JSON file "result.json"
my $jsonResult = "result.json"; # JSON result to get all info concerning CRISPR-Cas and sequences JSON files
open (JSONRES, ">$jsonResult") or die "open : $!";
print JSONRES "{\n";

my $dateJSON = $epoch; #$start_day."/".$start_month."/".$start_year."_".$start_hour.":".$start_min.":".$start_sec;
my $cmdLine = "perl $0 @ARGV";
my $jsonLineRes = "\"Date\":\"".$dateJSON."\",\n\"Version\":\"".$version."\",\n\"Command\":\"".$cmdLine."\",\n\"Sequences\":\n[";
##

if($logOption){
	print LOG "$cmdLine\n---> Results will be stored in $ResultDirFinal\n\n"; # the main command line and results path
	#[$start_hour:$start_min:$start_sec] 
}
if($quiet){
}
else{
	print " ---> Results will be stored in $ResultDirFinal\n\n"; #[$start_hour:$start_min:$start_sec]
}

#my $totalNumberOfCrisprs=0; #LK
my $allFoundCrisprs = 0; #DC
my $allCrisprs = 0;
my $nbrAllCas = 0; # Total nbrCas
my $currentRepository = getcwd();

#actualMetaOptionValue
my $actualMetaOptionValue = $metagenome;

#my header CCS to write title of CSS tsv file
my $resultsCRISPRCasSummary1 = $ResultDir."/CRISPR-Cas_summary.tsv";
my $headerCCS = "Sequence(s)\tCRISPR array(s)\tNb CRISPRs\tEvidence-levels\tCas cluster(s)\tNb Cas\tCas Types/Subtypes\n";
open (RESULTCCSONE, ">$resultsCRISPRCasSummary1") or die "open : $!";  #Open $resultsCRISPRCasSummary1
print RESULTCCSONE $headerCCS;
close(RESULTCCSONE);

#my Header for CRISPR-Cas_clusters.tsv
if($clusteringThreshold and $launchCasFinder){
  my $clusterResultsHead = $ResultDir."/CRISPR-Cas_clusters.tsv";
  my $headerCCC = "Sequence\tCluster Id\tNb CRISPRs\tNb Cas\tStart\tEnd\tDescription\n";
  open (RESULTCCC, ">$clusterResultsHead") or die "open : $!"; 
  print RESULTCCC $headerCCC;
  close (RESULTCCC);
 
}

#my Header for smallArraysReclassification.xls
if ($classifySmall){
  my $smallArrays = $ResultDir."/smallArraysReclassification.xls";
  open (SMALL, ">$smallArrays") or die "open : $!"; 
  print SMALL "CRISPR_Id\tConsensus_Repeat\tFormer_Evidence-level\tNew_Evidence-level\n";
  close (SMALL);
}

#sequence version
my $sequenceVersion = 0;

while($seq = $seqIO->next_seq()){  # DC - replace 'next_seq' by 'next_seq()'
  
  $inputfileCount++;
  $inputfile = $basename;

  my $start_run = time();
  # DC - 05/2017 - print $inputfile and seq->id
  my $seqID1 = $seq->id; # ID of sequence
  my $seqDesc = $seq->desc; # Description of sequence
  if($seqDesc eq ""){ $seqDesc = "Unknown"; }# If description is missing, the sequence will be labeled as "Unknown"
  
  $sequenceVersion = $seqID1;
  my $seqLength = $seq->length(); #Sequence size

  # -meta option as given by user
  $metagenome=$actualMetaOptionValue;

  # DC - when $seqLength < 100000 bases, consider it as a metagenomic dataset for Prodigal (if option -meta was not set)
  if(!$metagenome and ($seqLength < 100000) ){
  	$metagenome=1;
  }
  #else{
  #	$metagenome=$actualMetaOptionValue;
  #}
  

  #NV condition for mss option
  if ($seqLength >= $seqMinSize){
  	my $globalSeq = $seq->seq; # whole sequence to calculate global AT% # DC modif to perform better

  	my $globalAT = atpercent($globalSeq); # calculate AT% of running sequence # DC modif to perform better


  	my @seqID1 = (); # DC
  	if ($seqID1 =~ /\|/){
    		@seqID1 = split (/\|/,$seqID1);
    		$seqID1 = pop(@seqID1);
    		@seqID1 = ();
    		@seqID1 = split (/\./,$seqID1);
    		$seqID1 = $seqID1[0];
  	}
  	else{
    		if($seqID1 =~ /\./){
      			@seqID1 = split (/\./,$seqID1);
      			$seqID1 = $seqID1[0];
    		}
		
  	}
	# Add sequence number to variable $seqID1
	#$seqID1 = $seqID1."_".$inputfileCount;
  
  	# DC - 05/2017 - change value of $inputfile
  	$inputfile = $seqID1.".fna";
  
  	# DC - 05/2017 - To change the way of write sequence 
  	my $seqO = Bio::SeqIO->new(-file=>">$inputfile", -format=>'Fasta'); # DC
  	$seqO->write_seq($seq); # DC
  
  	my $nbrcris = 0;
  	my $OneSpacerCris_nbr = 0;  

  	# execute the mkvtree function to indexate the input sequence file
  	my $indexname = $inputfile;
  	my @indexname = split (/\./,$indexname);
  	print STDERR "\nSequence number $inputfileCount..";
  	#print STDERR "\n*** your results files will be in the $outdir-$inputfileCount directory ***\n"; 
  
  	callmkvtree($inputfile,$indexname);
  	# DC - 05/2017 print 
  	if($quiet){}
  	else{
  		print "  ( Input file: $inputfile, Sequence ID: $seqID1, Sequence name = $seqDesc )\n"; # DC
  	}
  	# execute the first vmatch search
  
  	# Modification DC - 05/05/2017
  	my @vmatchoptions = ""; #qw(-l 23 25 60 -e 1 -s leftseq -evalue 1 -absolute -nodist -noevalue -noscore -noidentity -sort ia -best 1000000 -selfun ./sel392v2.so 55); #DC
  	#my @vmatchoptions = qw(-l 23 25 60 -e 1 -s leftseq -evalue 1 -absolute -nodist -noevalue -noscore -noidentity -sort ia -best 1000000 -selfun sel392v2.so 55);
  
  	#my $currentRepositoryForSo = getcwd(); #get current repository to use sel392v2.so
  	#my ($vmatch_hour,$vmatch_min,$vmatch_sec) = Now();
  	# DC - set Vmatch options in function of parameters   ### replace -sort ia by -sort ida
  	#print " $repeatsQuery exists ?\n";
  	if(-e $so){
	  
	  if(!$mismOne){
	     if(-e $repeatsQuery){
		#print " $repeatsQuery exists !!!!\n";
		@vmatchoptions = ("-l", $M1, "-s", "leftseq", "-evalue", "1", "-absolute", "-nodist","-noevalue", "-noscore", "-noidentity", "-sort", "ia", "-q", "$repeatsQuery", "-best", 1000000, "-selfun", "$so", $M2);
			
		}
	     else{
		@vmatchoptions = ("-l", $M1, $S1, $S2,"-s", "leftseq", "-evalue", "1", "-absolute", "-nodist","-noevalue", "-noscore", "-noidentity", "-sort", "ia", "-best", 1000000, "-selfun", "$so", $M2);
	     }
	  }else{
	     if(-e $repeatsQuery){
		#print " $repeatsQuery exists !!!!\n";
	  	@vmatchoptions = ("-l", $M1, "-e $mismOne -s", "leftseq", "-evalue", "1", "-absolute", "-nodist","-noevalue", "-noscore", "-noidentity", "-sort", "ia", "-q", "$repeatsQuery", "-best", 1000000, "-selfun", "$so", $M2);
	     }
      	     else{
		@vmatchoptions = ("-l", $M1, $S1, $S2,"-e $mismOne -s", "leftseq", "-evalue", "1", "-absolute", "-nodist","-noevalue", "-noscore", "-noidentity", "-sort", "ia", "-best", 1000000, "-selfun", "$so", $M2);
	     }
	  }
  	}
  	else{
  		print "The shared object file ($so) must be available in your current directory. Otherwise, you must use option -soFile (or -so)! \n\n";
	
		if($logOption){
			print LOG "The shared object file ($so) must be available in your current directory. Otherwise, you must use option -soFile (or -so)! \n\n";  #[$vmatch_hour:$vmatch_min:$vmatch_sec] 
		}
        	exit(1);
  	}

  	push(@vmatchoptions,$indexname); # DC - replace 


  	push(@vmatchoptions, " > vmatch_result.txt");

  	if($logOption){
  		print LOG "\n vmatch @vmatchoptions\n"; # print in Logfile DC replaced vmatch by vmatch2 , [$vmatch_hour:$vmatch_min:$vmatch_sec] 
  	}
  	#Modification DC - 05/05/2017
  	#makesystemcall("./vmatch " . join(' ',@vmatchoptions)); #DC
  	makesystemcall("vmatch " . join(' ',@vmatchoptions)); #LK - DC replaced vmatch by vmatch2
  	#makesystemcall("cp vmatch_result.txt CHECKvmatch_result.txt");
  
  	my @rep = trans_data("vmatch_result.txt");
  

  	# fill in tabseq table : table containing the begin and end positions 
  	# of sequences candidates ton contain a crispr
  	# my @temp = split(/\./, $inputfile);
  	my $RefSeq = $seqID1; #$temp[0];
  	if($#rep >=0)
  	{
     		($nbrcris, $OneSpacerCris_nbr)=  write_clusters($RefSeq,@rep);
     

     		create_recap($RefSeq, $nbrcris, $OneSpacerCris_nbr,$ResultDir);
     		#$totalNumberOfCrisprs+=$OneSpacerCris_nbr; #LK
     		# DC
    		#$totalNumberOfCrisprs = $nbrcris;
     
  	}
  	else{ create_recap($RefSeq, 0, 0,$ResultDir); }
  	#$totalNumberOfCrisprs+=$nbrcris; #LK

  	#before changing directory/path
  	my $actual_pathHome = getcwd();
  	if($logOption){
		print LOG "Actual path in Home is: $actual_pathHome\n";
		#print "Actual path in Home is: $actual_pathHome\n";
  	}  

  	chdir ".."; #LK  #removed by DC
  	unlink <..\/crispr_result_*>;
  	unlink <..\/seq_v*>;
  	unlink <..\/$inputfile.*> ; 
  	unlink <..\/spacersseq_v*>;
    unlink <..\/fastaMuscle_spacersseq_v*>;

  	unlink '../stdout';
  	unlink '../vmatch_result.txt';
  	unlink '../alignDR_Spacer.needle';
  	unlink '../DR';
  	unlink <*_CRISPRs> ;
  	#unlink "../$inputfile"; # DC - 12/05/2017 to keep fasta file until gff and json files are done
  	chdir "..";
  	#unlink '*.index';
	
	makesystemcall("rm -f $inputfile.index");

  	## DC - move makeGFF and make JSON functions
  
  	# LK - DC
  	my ($gffFilename,@idDir)=makeGff($ResultDir,$inputfile,$seqDesc,$nbrcris,$OneSpacerCris_nbr);

  	# DC - 05/2017 - Make Json file thanks to GFF
  	#print "$gffFilename\n"; # DC
  	my $jsonFile = makeJson($gffFilename,$ResultDir,$RefSeq);
 

  	# DC - 05/2017 - CRISPRs analysis
  	if($quiet){}
  	else{
  		print "Nb of CRISPRs in this sequence = $nbrcris\n";
  	}
  	if($logOption){
    		print LOG "Nb of CRISPRs in this sequence $RefSeq = $nbrcris\n";
  	}

  	my ($analysis, $foundCrisprs) = crisprAnalysis($ResultDir,$RefSeq,$nbrcris,$seqDesc,@idDir);
  
  	$allFoundCrisprs += $foundCrisprs; #increment for each found CRISPR in CRISPRdb
  	$allCrisprs += $nbrcris; #all Crisprs
    
  	#report end time for CRISPRs analysis
  	my $end_runCRISPR = time();
  	my $runtimeCRISPR = $end_runCRISPR - $start_run;
    
  	# DC - 06/2017 - Integration of Cas analysis
  	my ($casfile, $jsonCAS, @tabCRISPRCasClustersA);
  	$casfile = "";
  	my $nbrCas = 0;

  	if($rcfowce){
	if($launchCasFinder and ($nbrcris>0)){
		($nbrCas, $casfile, $jsonCAS, @tabCRISPRCasClustersA) = casFinder($ResultDir,$inputfile,$seqDesc,$RefSeq,$nbrcris,$kingdom);
	}

  	}
  	elsif($launchCasFinder and ($rcfowce == 0)){
		($nbrCas, $casfile, $jsonCAS, @tabCRISPRCasClustersA) = casFinder($ResultDir,$inputfile,$seqDesc,$RefSeq,$nbrcris,$kingdom);
  	}
  	#print "NBR CAS = $nbrCas\n";
  	if(! $nbrCas){ $nbrCas = 0; }

  	if($quiet){}
  	else{
  		print "Nb of Cas in this sequence = $nbrCas\n";
  	}
  	if($logOption){
    		print LOG "Nb of Cas in this sequence ($RefSeq) = $nbrCas\n";
  	}
  
  	$nbrAllCas += $nbrCas; # all Cas

  	#if ($launchCasFinder){
  	#  $casfile = casFinder($ResultDir,$inputfile,$seqDesc,$RefSeq,$nbrcris,$kingdom,$casfinder);

  	#  print "Cas file = $casfile\n";

  	#}

  	#HTML page
  	if($html){
    		my $htmlPage = "";
    		$htmlPage = makeHtml($gffFilename,$casfile,$ResultDir,$RefSeq,$seqDesc,$seqLength,$globalAT,$nbrcris,$OneSpacerCris_nbr);
    		if($quiet){}
    		else{ print "HTML page ($htmlPage) created/updated\n"; }
  	}

  	#create a directory GFF and move all GFFs to this directory
  	#my $directoryGFFfiles = $ResultDir."/GFF";  # Directory for GFFs
  	#mkdir $directoryGFFfiles unless -d $directoryGFFfiles;
	if(-d $ResultDir."/GFF"){ makesystemcall("cp $gffFilename $ResultDir/GFF "); }

  	#Move CRISPRproperties, DRs and Spacers to specific directories
  	my $input_spacers = $ResultDir."/".$RefSeq."/Spacers_"; # Path to the Spacers files
  	my $input_drs = $ResultDir."/".$RefSeq."/DRs_"; # Path to the DRs files
  	my $input_properties = $ResultDir."/".$RefSeq."/".$RefSeq; # Path to the properties files

  	my $directorySpacers = $ResultDir."/".$RefSeq."/Spacers";  # Directory for spacers
  	my $directoryDRs = $ResultDir."/".$RefSeq."/DRs";  # Directory for DRs
  	my $directoryCRISPRproperties = $ResultDir."/".$RefSeq."/CRISPRproperties";  # Directory for CRISPRproperties
  	mkdir $directorySpacers unless -d $directorySpacers;
  	mkdir $directoryDRs unless -d $directoryDRs;
  	mkdir $directoryCRISPRproperties unless -d $directoryCRISPRproperties;
  
  	if ( (-e $input_spacers."1") and (-e $input_drs."1") ){
    		makesystemcall("cp -r $input_spacers* $directorySpacers");
    		makesystemcall("cp -r $input_drs* $directoryDRs");
    		makesystemcall("cp -r $input_properties* $directoryCRISPRproperties");
 	}

  	if(-e $input_properties."_CRISPRs"){
    		my $reportProperties = $input_properties."_CRISPRs";
    		makesystemcall("cp -r $reportProperties $directoryCRISPRproperties");
  	}

  	my $directoryProperties = $ResultDir."/".$RefSeq;
	#my $newProperties = $directoryProperties."_properties"; #before NV
	my $newProperties = ""; #NV
	if (-e $ResultDir."/CRISPRFinderProperties" and -d $ResultDir."/CRISPRFinderProperties") {
		$newProperties = $ResultDir."/CRISPRFinderProperties/".$RefSeq."_properties";
	} #NV
  	if(-d $directoryProperties){
    		makesystemcall("cp -r $directoryProperties $newProperties"); # rename properties for each sequence
  	}

  	#create a directory Properties and move all sequenceProperties directories in one
  	#my $directoryPropertiesFinal = $ResultDir."/CRISPRFinderProperties";  # Directory for Properties
  	#mkdir $directoryPropertiesFinal unless -d $directoryPropertiesFinal;
  	#makesystemcall("cp -r $newProperties $directoryPropertiesFinal");

  	# End Move properties, DRs and Spacers
  

  	# DC - 06/2017 - write sequence information
  	my $end_run = time();
  	my $runtime = $end_run - $start_run;
    	my $runtimeCas = $runtime - $runtimeCRISPR;#get calculation time for Cas analysis

  	if($logOption){
  		print LOGSEQ "$seqID1\t$seqDesc\t$seqLength\t$globalAT\t$nbrcris\t$runtime\t$runtimeCRISPR\t$runtimeCas\t$statusLOG[7]\n"; #some statistics
  	}

  	#$jsonLineSeq .= "{\n\"Id\":\"".$RefSeq."\",\n\"Description\":\"".$seqDesc."\",\n\"Length\":".$seqLength."\n},\n"; #data to be written in JSONSEQ file

  	##JSON result
  	my ($catJSONcrispr, $catJSONcas);

  	#my $jsonFileCrispr = $ResultDir."/".$RefSeq."_crispr.json";
  	if($jsonFile){
  		$catJSONcrispr = $jsonFile;
  	}
  	else{
		$catJSONcrispr = "[]";
  	}

  	#my $jsonFileCas = $ResultDir."/".$RefSeq."_cas.json";
  	if($jsonCAS){
  		$catJSONcas = $jsonCAS;
  	}
  	else{
		$catJSONcas = "[]";
  	}

  	if($nbrcris == 0){
		$catJSONcrispr = "[]";
  	}
  	if ( ! $jsonCAS){
  		$catJSONcas = "[]";
  	} 
	#get information from @tabCRISPRCasClustersA
  	my $tempSumDataCcC = "";
  	if ($launchCasFinder){
		$tempSumDataCcC = "@tabCRISPRCasClustersA";
  	}
  	else{
		$tempSumDataCcC = "NA";
  	}
#$sequenceVersion
  	$jsonLineRes .= "{\n\"Id\":\"".$RefSeq."\",\n\"Version\":\"".$sequenceVersion."\",\n\"Description\":\"".$seqDesc."\",\n\"AT\":".$globalAT.",\n\"Length\":".$seqLength.",\n\"Summary_CRISPR-Cas\":\"".$tempSumDataCcC."\",\n\"Crisprs\":".$catJSONcrispr.",\n\"Cas\":".$catJSONcas."\n},\n";


  	#put analyzed sequence in a folder named: "analyzedSequences"
  	my $analyzedSequences = $ResultDir."/analyzedSequences";  # Directory for FASTAs
  	mkdir $analyzedSequences unless -d $analyzedSequences;
  	if(-e $inputfile){
    		makesystemcall("cp -r $inputfile $analyzedSequences");
  	}
 } #NV end condition for mss option  
} # end while (my $seq = $seqIO->next_seq)

# End and close file JSONRES
$jsonLineRes .= "]\n";
$jsonLineRes =~ s/,\n]/\n]/; # change end of file to better fit JSON format
$jsonLineRes .= "}\n";

#$jsonLineRes =~ s/,/,\n/; # make file more readable

print JSONRES $jsonLineRes;
close (JSONRES);
# move JSONRES in ResultDir
if(-e $jsonResult){
  makesystemcall("cp -r $jsonResult $ResultDir");
}

#create a directory JSONfiles and move all JSONs to this directory
#my $directoryJSONfiles = $ResultDir."/JSON_files";  # Directory for JSONs
#mkdir $directoryJSONfiles unless -d $directoryJSONfiles;
#if(-e $ResultDir."/result.json"){
#  makesystemcall("cp -r $ResultDir/*.json $directoryJSONfiles");
#}

#fullReport
if ($launchCasFinder){
	if($writeFullReport){
	  my $fullRep = fullReport($ResultDir);
	  #print "CRISPR-CAS vicinity Report = $fullRep\n";
	}
}

#OrientationCount
if($dirRepeat){
    my $orientationCountFile = countOrientation($ResultDir,$allCrisprs);
    if($quiet){}
    else{ print "Orientations count file created: $orientationCountFile\n\n"; }
}

#copy files
my $resCCSummary = $ResultDir."/CRISPR-Cas_summary.tsv";
my $resCCSummaryX = $ResultDir."/CRISPR-Cas_summary.xls";
if(-e $resCCSummary){ makesystemcall("cp $resCCSummary $resCCSummaryX"); }

my $resCRISPRs = $ResultDir."/Crisprs_REPORT.tsv";
my $resCRISPRsX = $ResultDir."/Crisprs_REPORT.xls";
if(-e $resCRISPRs){ makesystemcall("cp $resCRISPRs $resCRISPRsX"); }

my $resCas = $ResultDir."/Cas_REPORT.tsv";
my $resCasX = $ResultDir."/Cas_REPORT.xls";
if(-e $resCas){ makesystemcall("cp $resCas $resCasX"); }

my $clustersResults = $ResultDir."/CRISPR-Cas_clusters.tsv";
my $clustersResultsX = $ResultDir."/CRISPR-Cas_clusters.xls";
if(-e $clustersResults){ makesystemcall("cp $clustersResults $clustersResultsX"); }

#create a directory TSVfiles and move all TSVs to this directory
my $directoryTSVfiles = $ResultDir."/TSV";  # Directory for TSVs
mkdir $directoryTSVfiles unless -d $directoryTSVfiles;
if(-e "$ResultDir/Crisprs_REPORT.tsv"){
  makesystemcall("cp $ResultDir/*.tsv $directoryTSVfiles");
}
if(-e "$ResultDir/Crisprs_REPORT.xls"){
  makesystemcall("cp $ResultDir/*.xls $directoryTSVfiles");
}

#create a directory FASTAfiles and move all rawFASTAs to this directory
my $directoryFASTAfiles = $ResultDir."/rawFASTA";  # Directory for FASTAs
mkdir $directoryFASTAfiles unless -d $directoryFASTAfiles;
#if(-e "$ResultDir/rawCRISPRs.fna"){
  #makesystemcall("mv $ResultDir/rawCRISPRs.fna $directoryFASTAfiles");
#}
#if(-e "$ResultDir/rawCas.fna"){
  #makesystemcall("mv $ResultDir/rawCas.fna $directoryFASTAfiles");
#}

my $analyzedSequences2 = $ResultDir."/analyzedSequences";  # Directory for analyzed FASTAs
if(-d $analyzedSequences2){
  makesystemcall("cp -r $analyzedSequences2 $directoryFASTAfiles");
}

#create a directory Properties and move all sequenceProperties directories in one
my $directoryProperties = $ResultDir."/CRISPRFinderProperties";  # Directory for Properties
#mkdir $directoryProperties unless -d $directoryProperties; #before NV
#if(-e "$ResultDir/*_properties"){
  #makesystemcall("cp -r $ResultDir/*_properties $directoryProperties"); #before NV
#}

#create a directory GFF and move all GFFs to this directory
#my $directoryGFFfiles = $ResultDir."/GFF";  # Directory for GFFs # DC
#mkdir $directoryGFFfiles unless -d $directoryGFFfiles; # DC
#if(-e "$ResultDir/*.gff"){
  #makesystemcall("cp -r $ResultDir/*.gff $directoryGFFfiles"); #before NV
  #makesystemcall("find $ResultDir -name *.gff -exec mv {} $directoryGFFfiles \\;"); #NV DC
#}

## Option Keep to remove (by default) all secondary files or keep them
if($keep){
  if($useProkka){
  	if(! $quiet){print "Secondary folders/files (Prokka, CasFinder, rawFASTA, CRISPRFinderProperties) have been created\n";}
  }
  else{
	if(! $quiet){print "Secondary folders/files (Prodigal, CasFinder, rawFASTA, CRISPRFinderProperties) have been created\n";}
  }
}
else{
  if(! $quiet){print "Deleting secondary folders/files...\n";}
  if($launchCasFinder){
    if($useProkka){
    	makesystemcall("rm -rf $ResultDir/Prokka") if (-d "$ResultDir/Prokka"); #rmdir $ResultDir."/Prokka";
    }
    else{ 
    	makesystemcall("rm -rf $ResultDir/Prodigal") if (-d "$ResultDir/Prodigal"); #rmdir $ResultDir."/Prodigal";
    }

    makesystemcall("rm -rf $ResultDir/CasFinder") if (-d "$ResultDir/CasFinder"); #rmdir $ResultDir."/CasFinder"; #
  }

  makesystemcall("rm -rf $ResultDir/CRISPRFinderProperties") if (-d "$ResultDir/CRISPRFinderProperties"); # rmdir $ResultDir."/CRISPRFinderProperties"; #
  makesystemcall("rm -rf $ResultDir/rawFASTA") if (-d "$ResultDir/rawFASTA"); #rmdir $ResultDir."/rawFASTA"; #
}

if($crisprdb eq ""){
  print "\nAll CRISPRs = $allCrisprs\n";
}
elsif(-e $crisprdb){
  print "\nNb CRISPRs matching CRISPRdb / All CRISPRs = $allFoundCrisprs / $allCrisprs\n";
}

#print Cas
if($launchCasFinder){
  print "All Cas = $nbrAllCas\n";
}

## Time - end
#my ($end_hour,$end_min,$end_sec) = Now();
#my($end_year,$end_month,$end_day, $end_hour,$end_min,$end_sec) = Today_and_Now();
##

#my ($D_y,$D_m,$D_d, $Dh,$Dm,$Ds) =
 #     Delta_YMDHMS($start_year,$start_month,$start_day, $start_hour, $start_min, $start_sec,
 #                  $end_year, $end_month, $end_day, $end_hour,$end_min,$end_sec);
 
print "\nThank you for using $0! Thank you for your patience!\n"; # [$end_hour:$end_min:$end_sec] 
#print "\n[$end_hour:$end_min:$end_sec] The script lasted: ".$D_y." year(s) ".$D_m." month(s) ".$D_d." day(s) , ".$Dh." hour(s) ".$Dm." minute(s) ".$Ds." second(s)\n";
my $endEpoch = time(); #my $epoch = time();
my $minutesEpoch = ($endEpoch - $epoch) / 60;

print "\nThe script lasted: ".($minutesEpoch/60)." hour(s) ".$minutesEpoch." minute(s) ".($endEpoch - $epoch)." second(s)\n";


if($logOption){
  #print LOG "\n[$end_hour:$end_min:$end_sec] The script lasted: ".$D_y." year(s) ".$D_m." month(s) ".$D_d." day(s) , ".$Dh." hour(s) ".$Dm." minute(s) ".$Ds." second(s)\n";

  close (LOG);
  close (LOGSEQ);
}
#close HTML page if exists
if($html){
  print HTML "</div>";
  print HTML "</body>";
  print HTML "</html>";

  close (HTML);
}

# Move LOG files in LOGs directory

if($logOption){
  my $directoryLog = $ResultDir."/LOGs"; # directory  to store LOG files
  unless(-d $directoryLog){ mkdir $directoryLog or die "$0: I can not create the folder $directoryLog: $!\n" }

  if ( (-e $logfile) and (-e $logSeq) ){
    makesystemcall("cp $logfile $directoryLog");
    makesystemcall("cp $logSeq $directoryLog");
  }
}

# Move index.html and copy crispr.css into Visualization directory
if($html){
  my $directoryViz = $ResultDir."/Visualization"; # directory  to store basic Visualization files
  unless(-d $directoryViz){ mkdir $directoryViz or die "$0: I can not create the folder $directoryViz: $!\n" }

  makesystemcall("cp $htmlFile $directoryViz");

  if(-e "crispr.css"){
  	makesystemcall("cp crispr.css $directoryViz");
  }
  elsif(-e "supplementary_files/crispr.css"){
	makesystemcall("cp supplementary_files/crispr.css $directoryViz");
  }
  elsif(-e $cssFile){
	makesystemcall("cp $cssFile $directoryViz/crispr.css");
  }
}

## Deleting temporaryFolder and move to ResultDirFinal

  if(-d "$outdir-$inputfileCount"){  
  	`rm -Rf $outdir-$inputfileCount &>/dev/null`;
  }
  #`mv $ResultDir/$outdir $ResultDir/$outdir-$inputfileCount`; #LK - removed by DC

  ## DC - $outdir will be moved to $ResultDir after
  #rmdir "$outdir-$inputfileCount"; # Only if empty
  #`mv $ResultDir $ResultDirFinal`;#rmdir $outdir;
  #makesystemcall("mv $ResultDir $ResultDirFinal") if (-d $ResultDir);

# Remove '.index' files
makesystemcall("rm -f *.index");

# Remove macsyOutput if option -quiet is set
if($quiet){
	makesystemcall("rm -f macsyfinderOutput");
}

# // LK
# DC - End of main script

#TODO make a GBROWSE file too # Maybe in a next future.
# http://nbase.biology.gatech.edu/gbrowse2/general_help.html#upload

#------------------------------------------------------------------------------
################ /CRISPRCasFinder functions

#------------------------------------------------------------------------------
# DC - 07/2017 - makeHTML (creates HTML file from GFF file and Cas report)
sub makeHtml
{
  my ($gff,$casFile,$ResultDir,$refSeq,$seqDesc,$seqLength,$globalAT,$nbrcris,$OneSpacerCris_nbr) = @_; # parameters

  my $input_spacer = $ResultDir."/".$refSeq."/Spacers_"; # Path to the Spacers files
  my $input_dr = $ResultDir."/".$refSeq."/DRs_"; # Path to the DRs files

  my $nbConfirmed = $nbrcris - $OneSpacerCris_nbr; #get the number of confirmed CRISPRs arrays
  
  
  my ($crisprID,$name,$orientation,$start,$end,$DRconsensus,$DRlength,$nbSpacers,$sequence,$leader);

  my @colors = ("Bisque", "Blue", "Blueviolet","Brown","Burlywood","Cadetblue","Chartreuse","Chocolate","Coral","Cornflowerblue","Crimson","Cyan","Darkcyan","Darkgoldenrod","Darkgray","Darkgreen","Darkkhaki","Darkmagenta","Darkolivegreen",
"Darkorange","Darkorchid","Darkred","Darksalmon","Darkseagreen","Darkslategray",
"deeppink","Deepskyblue","Moccasin","Dodgerblue","Firebrick","Lemonchiffon","Forestgreen",
"Fuchsia","Gainsboro","Gold","Goldenrod","Gray","Green","Lightgoldenrodyellow","Greenyellow","Hotpink",
"Indianred","Indigo","Mediumblue","Khaki","Lavender","Lawngreen","Lightblue",
"Lightcoral","Lightcyan","Darkturquoise","Lightgreen","Lightgrey","Lightpink","Lightsalmon","Lime",
"Lightseagreen","Lightskyblue","Lightslategray","Lightsteelblue","Limegreen",
"Magenta","Maroon","Mediumorchid","Mediumpurple","Aquamarine","Mediumseagreen",
"Mediumslateblue","Mediumspringgreen","Mediumturquoise",,"Navajowhite","Mediumvioletred","Mistyrose","Navy","Olive","Aliceblue","Olivedrab","Orange","Orangered","Orchid",
"Palegoldenrod","Palegreen","Paleturquoise","Palevioletred","Peachpuff","Peru","Pink","Plum",
"Powderblue","Purple","Red","Rosybrown","Royalblue","Salmon","Sandybrown","Seagreen",
"Sienna","Silver","Skyblue","Slateblue","Slategray","Springgreen","Steelblue","Tan","Teal","Thistle",
"Tomato","Turquoise","Wheat","YellowGreen");
  my $colorcount = 0;


  # Read and parse GFF file
  open(GFF,"<",$gff) or die("Could not open the GFF file $gff because $!\n");
  
  # Read and parse Cas report file
  if(-e $casFile){
    open(CAS,"<",$casFile) or die("Could not open the Cas report file $casFile because $!\n");
  }

  #my ($hour,$min,$sec) = Now();
  if($logOption){
  	#print LOG "\n[$hour:$min:$sec] Create HTML file (index.html) for visualization\n";
  }

  #open(HTML,">>",$htmlFile) or die("Could not open the file $htmlFile because $!\n");

  print HTML "<br/><br/><br/>";
  print HTML "<table class= 'primary_table' width = '100%'>";
  print HTML "<tr><th>Sequence: ", $refSeq,"</th></tr></table>\n";	 	
  print HTML "<br/><font COLOR= \#000099> Sequence description : </font>", $seqDesc;
  print HTML "<br/><font COLOR= \#000099> Length (bp): </font>", $seqLength,"<br/>\n";
  print HTML "<br/><font COLOR= \#000099> AT%: </font>", $globalAT,"<br/>\n";
  print HTML "<br/><font COLOR= \#000099> Number of CRISPRs candidates: </font>", $nbrcris,"<br/>\n";
  #print HTML "<br/><font COLOR= \#000099> - Number of Confirmed CRISPRs: </font>", $nbConfirmed,"<br/>\n";
  #print HTML "<br/><font COLOR= \#000099> - Number of Hypothetical CRISPRs: </font>", $OneSpacerCris_nbr,"<br/><br/><br/>\n";

  print HTML "<center><table class = 'sub' width = '80%'><tr><th>CRISPRs</th></tr></table></center><br/><br/>\n";

  my $leftFlank = "";
  my $rightFlank = "";

  my ($entropyDRs, $conservationSpacers, $meanSpacers, $ratioDRspacer,$leftAT,$rightAT);

  while(<GFF>) {
    chomp($_);

    if ($_ !~  m/^#/ && $_ =~  m/CRISPRCasFinder/) {
        
        my @tab = split(/\t/, $_);

	#### LeftFLANK (left flanking sequence)
	if($tab[2] eq "LeftFLANK"){
	  $start = $tab[3];
	  $end = $tab[4];
	  my @newTab = split(/;/, $tab[8]);
	  $sequence = $newTab[0];
	  $sequence =~ s/sequence=//;
  	  $leftAT = $newTab[1];
   	  $leftAT =~ s/at%=//;

  	  $leftFlank = "($flankingRegion bp [start=$start;end=$end], AT%=$leftAT): ".$sequence;
	  
	}

	#### main CRISPR
        if($tab[2] eq "CRISPR" || $tab[2] eq "PossibleCRISPR"){

	  my %hash = ();
	  $start = $tab[3];
	  $end = $tab[4];

	  $orientation = $tab[6];

	  my @newTab = split(/;/, $tab[8]);
	  for (my $i = 0; $i <= $#newTab; $i++) {
	    my @otherTab = split(/=/, $newTab[$i]);
	    $hash{$otherTab[0]} = $otherTab[1];
	  }
	  $DRconsensus = $hash{DR};
	  $DRlength = $hash{DR_length};
	  $nbSpacers = $hash{Number_of_spacers};
	  $crisprID = $hash{ID};
	  $name = $hash{name};

	  ## Get ID
	  my @tabIdCRISPR = split(/_/, $crisprID);
	  my $idNumber = pop(@tabIdCRISPR);


	  ## Entropy DRs
	  if(-e $input_dr.$idNumber){
	  	my $drFasta = fastaAlignmentMuscle($input_dr.$idNumber);
	  	$entropyDRs = entropy($drFasta);
	  }

	  ## Conservation Spacers
	  if(-e $input_spacer.$idNumber) {
	  	my $spacerFasta = fastaAlignmentMuscle($input_spacer.$idNumber);
		$conservationSpacers = sequenceAlignmentMuscle($spacerFasta);

		my @tabSpacerLength=();
		open (SPACER, "<".$input_spacer.$idNumber) or die "open : $!";  #Open Spacer file

		while (defined (my $lineSp = <SPACER>)){
			chomp($lineSp);
			my $lengthInfo = 0;
			if($lineSp !~ /^>/){
				#$lineSp =~ s/\s//g;   #Remove spaces if exist
				$lengthInfo = length($lineSp);
				push(@tabSpacerLength,$lengthInfo);
			}
		}
		close SPACER;

		##add Mean Spacers
		$meanSpacers = &average(\@tabSpacerLength);
		#$medianSpacers = median(@tabSpacerLength);
		##add ratio DR/Spacer
		$ratioDRspacer = $DRlength / $meanSpacers;

	  }

	  # Transform $crisprID
	  $crisprID =~ s/PossibleCrispr_//; #PossibleCrispr_  Crispr_
	  $crisprID =~ s/Crispr_//;

	  ## Evidence Level
	  my $eL = 0; # evidenceLevel (eL) = 4 means good CRISPR; eL = 3 means acceptable CRISPR; eL = 2 means bad CRISPR; eL = 1 means hypothetical CRISPR; 

	  if( ($entropyDRs >= 70) and ($conservationSpacers <= 8) and ($nbSpacers > 3) ){ $eL=4;}
	  if( ($entropyDRs >= 70) and ($conservationSpacers > 8) and ($nbSpacers > 3) ){ $eL=3;}
	  if( ($entropyDRs < 70) and ($nbSpacers > 3) ){ $eL=2;}
	  if( $nbSpacers <= 3 ){ $eL=1;}

	  #if( ($entropyDRs >= 70) and ($conservationSpacers <= 8) and ($nbSpacers > 3) and ($ratioDRspacer>=0.8) and ($ratioDRspacer<=1.2) ){ $eL=4;}
      	  #if( ($entropyDRs >= 70) and ($nbSpacers > 3) and (($conservationSpacers > 8) or ($ratioDRspacer<0.8) or ($ratioDRspacer>1.2) ) ){ $eL=3;}
      	  #if( ($entropyDRs < 70) and ($nbSpacers > 3) ){ $eL=2;}
      	  #if( $nbSpacers <= 3 ){ $eL=1;}

	  ## Repeat ID + CRISPRDirection
	  my @tabR = repeatIDandNb($DRconsensus);
      	  my $crisprDirection = repeatDirection($tabR[0]);
	  if($crisprDirection =~  m/^F/){
		$crisprDirection = "+";
	  }
	  elsif($crisprDirection =~  m/^R/){
		$crisprDirection = "-";
	  }
	  else{
		$crisprDirection = "ND";
	  }


	  print HTML "<h5>Candidate number $idNumber</h5>\n";
	  print HTML "<div class = box> \n";
	  print HTML "<B> <FONT COLOR= \#9900FF><a id =$crisprID> </a> CRISPR ID : </FONT> $crisprID </B> <UL>\n";
	  print HTML "<li><FONT COLOR= \#000099>Left flanking sequence: </FONT> $leftFlank </li>\n";
	  print HTML "\n <li><FONT COLOR= \#000099>CRISPR start position: </FONT> $start ----------";
	  print HTML "<FONT COLOR= \#000099> CRISPR end position: </FONT> $end ----------";
	  my $len = $end - $start;
	  print HTML "<FONT COLOR= \#000099> CRISPR length: </FONT> $len</li>\n";
	  print HTML "<li><FONT COLOR= \#000099> DR consensus: </FONT> $DRconsensus </li>\n<li><FONT COLOR= \#000099> DR length: </FONT> $DRlength";
	  print HTML "<FONT COLOR= \#000099> Number of spacers: </FONT> $nbSpacers </li>\n";
	  print HTML "<li><FONT COLOR= \#000099> Consensus Repeat ID (Nb in CRISPRdb): </FONT> $tabR[0] ($tabR[1]) </li>\n";
	  print HTML "<li><FONT COLOR= \#000099> CRISPRDirection: </FONT> $crisprDirection </li>\n";
	  print HTML "<li><FONT COLOR= \#000099> Potential orientation (AT%): </FONT> $orientation </li>\n";
  	  print HTML "<li><FONT COLOR= \#000099> Conservation of DRs based on Entropy: </FONT> $entropyDRs </li>\n";
	  print HTML "<li><FONT COLOR= \#000099> Conservation of Spacers: </FONT> $conservationSpacers </li>\n";
	  print HTML "<li><FONT COLOR= \#000099> Evidence Level: </FONT> $eL </li>\n";
	  print HTML "</UL>\n";
	  print HTML "<table class ='table2'  BGCOLOR=\#ffffff  width=70%>";

        }
	
	#### DR or Spacer
	if($tab[2] eq "CRISPRdr" || $tab[2] eq "CRISPRspacer"){
	  $start = $tab[3];
	  $end = $tab[4];
	  my @newTab = split(/;/, $tab[8]);
	  $sequence = $newTab[0];
	  $sequence =~ s/sequence=//;
	  my $seqType = "";

	  #repeats+spacers block
	  if($tab[2] eq "CRISPRdr"){
	    $seqType = "DR";
	    
	    #DC trying to color mismatches in red #COLORMISMATCH 
	    my @charSeqDRconsensus = split(//, $DRconsensus);
	    my @newDRpatternHtml = split(//, $sequence);
	    my $newChainDRhtml = "";

            print HTML "<tr> <td align = 'left'> <FONT size=2.5 color='\#000099'> $start </font></td>";

	    for(my $kh = 0; $kh < length($DRconsensus); $kh++){
	      if($charSeqDRconsensus[$kh] ne $newDRpatternHtml[$kh]){
		$newChainDRhtml .= "<FONT size=2.7 face='Courier New' color='\#FF0000'><b>$newDRpatternHtml[$kh]</b></FONT>";
	      }
	      else{
		$newChainDRhtml .= "<FONT size=2.7 face='Courier New'><b>$newDRpatternHtml[$kh]</b></FONT>";
   	      }
	    }
	    chomp($newChainDRhtml);
	    print HTML "<td align = 'left'><span style='background:Yellow;'>$newChainDRhtml</span></td>";

	    #print HTML "<td align = 'left'><span style='background:Yellow;'><FONT size=2.7 face='Courier New'> <b>$sequence</b> </FONT></span></td>";
	    
	  }
	  elsif($tab[2] eq "CRISPRspacer"){
	    $seqType = "Spacer";

	    my $col = $colors[$colorcount]; 
	    if($colorcount == $#colors){$colorcount = 0;}else{$colorcount++;}

	    print HTML "<td align = 'left'><span style='background:$col;'><FONT size=2.7 face= 'Courier New'>$sequence </FONT></span></td>";
	    print HTML "<td><FONT size=2.5 color = '\#000099'>$end</font></td></tr>";

	  }
	 
	}

	#### RightFLANK
	if($tab[2] eq "RightFLANK"){
	  $start = $tab[3];
	  $end = $tab[4];
	  my @newTab = split(/;/, $tab[8]);
	  $sequence = $newTab[0];
	  $sequence =~ s/sequence=//;
	  $rightAT = $newTab[1];
   	  $rightAT =~ s/at%=//;

	  $rightFlank = "($flankingRegion bp [start=$start;end=$end], AT%=$rightAT): ".$sequence;

	  if($newTab[2] eq "leader"){
	    $leader = 1;
	  }
	  else {
	    $leader = 0;
	  }

	  my $lastEnd = $start - 1;
	  print HTML "<td></td><td><FONT size=2.5 color = '\#000099'> $lastEnd </font></td></tr></table>\n";
	  print HTML "<ul><li><FONT COLOR= \#000099>Right flanking sequence: </FONT> $rightFlank </li></ul>\n";

	  print HTML "</div><br/><br/>\n";
	  

	}
	
    }

  }


  if(-e $casFile){
  
    print HTML "<center><table class = 'sub' width = '80%'><tr><th>Cas system(s) and genes</th></tr></table></center><br/><br/>\n";

    my $countSystem = 0;
    while(<CAS>) {
      chomp($_);

      if ($_ =~ m/####Summary system/ && $_ =~ m/$refSeq/) {
	 $countSystem++;
	 my $lineCas = $_;
	 $lineCas =~ s/####Summary system//; 

	 print HTML "<div class = box> \n";
	 print HTML "<ul><li><FONT COLOR= \#000099>System number $countSystem: </FONT> $lineCas </li></ul>\n";
      }

    }

    # end of div
    print HTML "</div>";
    print HTML "<br/><br/>";


  }

  close(GFF);

  if(-e $casFile){
    close(CAS);
  }


  return "index.html";
}

#------------------------------------------------------------------------------
# DC
#############################################################################
# Function allowing to parse MacSyFinder JSON file (with module JSON::Parse)
# and to parse GFF3 file generated with Prodigal or Prokka on given genomes.
# Results: list of Cas systems and Cas genes corresponding to analyzed genome
# with positions and orientations of genes, information relative to Cas,
# and other information (annotations)
# October 2022: changes regarding MacSyFinder v2 and CASFinder v3
#############################################################################
sub casFinder
{
  my($ResultDir,$inputfile,$seqDesc,$RefSeq,$nbrcris,$kingdom)=@_;

  my $repProkka = "";
  if($useProkka){ $repProkka = $ResultDir."/prokka_".$RefSeq; }
  else{ $repProkka = $ResultDir."/prodigal_".$RefSeq; }

  my $casDir =  $ResultDir."/casfinder_".$RefSeq."/"; 
 
  my $jsonCAS = ""; # variable allowing to write JSON file related to Cas

  my $nbCas = 0; # nb Cas for each search
  my $casdb = "";
  my $default = 0;
  my $addToMacSy = "";
  #manage $definition info
  # if( ($definition eq "General") or ($definition eq "general") or ($definition eq "G") or ($definition eq "g") ){
	# $casdb = $casfinder."/DEF-Class-2.0.3/";  # DEF-General-2.0 replaced by DEF-Class-2.0.2 then replaced by DEF-Class-2.0.3
	# $addToMacSy .= "--min-genes-required General-Class1 1 "; # General-1CAS replaced by General-Class1
	# $addToMacSy .= "--min-genes-required General-Class2 1 "; # General-3CAS replaced by General-Class2
	# #$default = 1;
  # }
  # elsif( ($definition eq "Typing") or ($definition eq "typing") or ($definition eq "T") or ($definition eq "t") ){
  	# $casdb = $casfinder."/DEF-Typing-2.0.3/"; # DEF-Typing-2.0 replaced by DEF-Typing-2.0.2 then replaced by DEF-Typing-2.0.3
	# #$addToMacSy .= "--min-genes-required CAS 1 --min-genes-required CAS-TypeI 1 --min-genes-required CAS-TypeII 1 --min-genes-required CAS-TypeIII 1 --min-genes-required CAS-TypeIV 1 --min-genes-required CAS-TypeV 1 ";
	# #$addToMacSy .= "--min-genes-required CAS-TypeVI 1 ";
  # }
  # elsif( ($definition eq "SubTyping") or ($definition eq "subtyping") or ($definition eq "S") or ($definition eq "s") ){
  	# $casdb = $casfinder."/DEF-SubTyping-2.0.3/"; # DEF-SubTyping-2.0 replaced by DEF-SubTyping-2.0.2 then replaced by DEF-SubTyping-2.0.3
  # }
  # else{
  	# #$casdb = $casfinder."/DEF-General-2.0/";
	# #$default = 1;
  # }
  
  # $casdb =~ s/\/\//\//; #DC - replace '//' by '/'

  # my $profiles = $casfinder."/CASprofiles-2.0.3/"; # CASprofiles-2.0 replaced by CASprofiles-2.0.2 then replaced by CASprofiles-2.0.3
  # $profiles =~ s/\/\//\//; #DC - replace '//' by '/'

  my $results = $ResultDir."/Cas_REPORT.tsv"; # former name: $ResultDir."/Cas_".$RefSeq.".txt";
  my $allCas = $ResultDir."/CRISPR-Cas_systems_vicinity.tsv"; # to store all individual CRISPR-Cas systems

  #my ($hour,$min,$sec) = Now(); 

  #open CCS.tsv file
  my $resultsCRISPRCasSummary = $ResultDir."/CRISPR-Cas_summary.tsv";
  open (RESULTCCS, ">>$resultsCRISPRCasSummary") or die "open : $!";  #Open $resultsCRISPRCasSummary
  my @tabCCScas = (); # table which will contain all Cas systems of the given sequence
  my %hashCountCas = (); # hash table to count Cas systems

  ####### work on CRISPR-Cas clusters and tracrRNA -Begin

	  my @tabCRISPRCasClusters = (); # table which will contain all CRISPR-Cas clusters according to a given threshold
	  my %hashCasBegin = ();
	  my %hashCasEnd = ();
	  my %hashCrisprBegin = ();
	  my %hashCrisprEnd = ();

	  my %hashCrisprDRcons = (); # to store the consensus DR corresponding to each CRISPR array
	  my $limitSide = 3000; # limit allowed in both side of Cas system or neighboring CRISPR
	  my $matchingCasType = "CAS-TypeII";

	  my @otherTabCas = ();
	  my @otherTabCrispr = ();

	  #my $clusterCrisprCas = "";
	  my $clusterResults = $ResultDir."/CRISPR-Cas_clusters.tsv";

  ####### work on CRISPR-Cas clusters and tracrRNA -End

  eval{
    #use JSON::Parse 'json_file_to_perl';
    #use JSON::Parse 'valid_json';
    
    my $prokkaProg = isProgInstalled("prokka");
    my $prodigalProg = isProgInstalled("prodigal");

    if($useProkka){ # if $useProkka
    if($prokkaProg){
      if(! $quiet) { print "\n\nprokka installation is.............OK \n"; }
    }
    elsif( (-e $userFAA) and (-e $userGFF) ){
      print "\n\nprokka format is used for GFF:$userGFF and proteome:$userFAA \n";
    }
    else
    {
      print "\n\n";
      print " _\n";
      print "/!\\ prokka is not installed ........ \n";
      print "___\n";
      print "Please install it by following the documentation provided here: https://github.com/tseemann/prokka  OR  http://www.vicbioinformatics.com/software.prokka.shtml\n";
      
      print "Otherwise, please retry without Cas option (-cas)\n\n";

      if($logOption){
        print LOG "\n/!\\ prokka is not installed ........ \n"; #[$hour:$min:$sec] 
        print LOG "\nPlease install it by following the documentation provided here: https://github.com/tseemann/prokka  OR  http://www.vicbioinformatics.com/software.prokka.shtml\n\n";
        print LOG "\nOtherwise, please retry without Cas option (-cas)\n\n";
      }

      #printhelpall($0);
      #exit 0;
      $launchCasFinder = 0;
      next;
    }
    }#end if $useProkka

    else{#$useProdigal
    if($prodigalProg){
      if(! $quiet) { print "\n\nprodigal installation is.............OK \n"; }
    }
    elsif( (-e $userFAA) and (-e $userGFF) ){
      print "\n\nprodigal format is used for GFF:$userGFF and proteome:$userFAA \n";
    }
    else
    {
      print "\n\n";
      print " _\n";
      print "/!\\ prodigal is not installed ........ \n";
      print "___\n";
      print "Please install it by following the documentation provided here: https://github.com/hyattpd/Prodigal\n";
      
      print "Otherwise, please retry without Cas option (-cas)\n\n";

      if($logOption){
        print LOG "\n/!\\ prodigal is not installed ........ \n"; #[$hour:$min:$sec] 
        print LOG "\nInstall it by following the documentation provided here: https://github.com/hyattpd/Prodigal\n\n";
        print LOG "\nOtherwise, please retry without Cas option (-cas)\n\n";
      }

      #printhelpall($0);
      #exit 0;
      $launchCasFinder = 0;
      next;
    }
    }

    my $macsyfinderProg = isProgInstalled("macsyfinder");
    if($macsyfinderProg){
      if (! $quiet) { print "macsyfinder installation is...........OK \n"; }
    }
    else
    {
      print "\n";
      print " _\n";
      print "/!\\ macsyfinder is not installed .........\n";
      print "___\n";
      print "Please install it by following the documentation provided here: https://github.com/gem-pasteur/macsyfinder\n";
      print "Otherwise, please retry without Cas option (-cas)\n\n";
      #printhelpall($0);
      #exit 0;
      $launchCasFinder = 0;
      next;
    }


    # Call Prokka
    #my ($prokka_hour,$prokka_min,$prokka_sec) = Now();
    my $options = "";
    if($quiet){ $options .= "--quiet "; }
    if($fast){ $options .= "--fast --rawproduct --norrna --notrna "; }
    if($metagenome){ $options .= "--metagenome "; }
    my $prokka = "prokka $options";
    $prokka .= "--cpus $cpuProkka --kingdom $kingdom --gcode $genCode --outdir $repProkka --prefix $RefSeq $inputfile"; # Prokka
    #print "$prokka\n";

    # Call Prodigal instead of Prokka (to be faster)
    # create repository if ($useProdigal)
    if($useProdigal and (! -e $userFAA) and (! -e $userGFF) ){
	unless(-d $repProkka){ mkdir $repProkka or die "$0: $repProkka can not be created: $!\n"; }
    }
    my $prodigal = "prodigal -i $inputfile -c -m -g $genCode -f gff -a $repProkka/$RefSeq.faa -o $repProkka/$RefSeq.gff "; 
    if($quiet){ $prodigal .= "-q "; }
    if($metagenome){ $prodigal .= "-p meta "; }
    else{ $prodigal .= "-p single "; }

    my $macsyfinder = "";
    my $proteome = ""; # variable to indicate if the user has entered a proteome or prokka generetated it

    if(-e $userFAA){
	$proteome = $userFAA;
    }
    else{
	$proteome = "$repProkka/$RefSeq.faa";
    }
	
    #$macsyfinder = "macsyfinder --sequence-db my_aa.fna --replicon-topology circular --db-type gembase --models CASFinder --accessory-weight 1 --exchangeable-weight 1 --coverage-profile 0.4 --redundancy-penalty 1 -w 4 -o $casDir ";
	$macsyfinder = "macsyfinder -w $cpuMacSyFinder --replicon-topology circular --models CASFinder --accessory-weight 1 --exchangeable-weight 1 --coverage-profile 0.4 --redundancy-penalty 1 --sequence-db $proteome --out-dir $casDir ";
	    
    #my ($macsyfinder_hour,$macsyfinder_min,$macsyfinder_sec) = Now();
    #if ( (-d $casfinder) and (-d $casdb) and (-d $profiles) ) {
	  #if($default){
	    #$macsyfinder = "macsyfinder -w $cpuMacSyFinder -d $casdb -p $profiles --sequence-db $proteome all --out-dir $casDir ";
	    #$macsyfinder = "macsyfinder -w $cpuMacSyFinder -d $casdb -p $profiles --sequence-db $proteome all --out-dir $casDir "; 
	    #$macsyfinder .= $addToMacSy; # --multi-loci General-CAS
	  #}
	  #else{
	    #$macsyfinder = "macsyfinder -w $cpuMacSyFinder -d $casdb -p $profiles --sequence-db $proteome all --out-dir $casDir "; #--multi-loci 'CAS*' --min-genes-required 'CAS*' 1";
	  #}
    #}
    #elsif($default){
        #$macsyfinder = "macsyfinder -w $cpuMacSyFinder --sequence-db $proteome all --out-dir $casDir"; #MacSyfinder without new definitions/profiles
    #}
    #else{  
    #	$macsyfinder = "macsyfinder -w $cpuMacSyFinder --sequence-db $proteome all --out-dir $casDir"; #MacSyfinder without new definitions/profiles
    #}
    #print "$macsyfinder\n";
    #print "METAGENOME = $metagenome ; ACTUAL VALUE = $actualMetaOptionValue\n\n\n";
    if($metagenome){ $macsyfinder .= " --db-type unordered"; } 
    elsif($actualMetaOptionValue == 0){ $macsyfinder .= " --db-type gembase"; } #--db-type ordered_replicon
    #elsif($gembase){ $macsyfinder .= " --db-type gembase"; }
    if($quiet){ 
	open (OUTMACSY, ">macsyfinderOutput") or die "open : $!";  # Open macsyfinderOutput (writing mode)
	$macsyfinder .= " >> macsyfinderOutput";
	if($logOption){
      		print LOG "\n MacSyFinder output when launched quieter\n"; #[$hour:$min:$sec]
        	open (OUTMACSY, "<macsyfinderOutput") or die "open : $!";  # Open macsyfinderOutput
		while (my $fileOutMacsy = <OUTMACSY>) {
         		print LOG "$fileOutMacsy";
    		} 
		print LOG "\n";
		close(OUTMACSY);
		system ("rm -f macsyfinderOutput");
	}
    } #&>> outputMacSyFinder.out   or -v

    my $tsv = $casDir."best_solution.tsv"; # the JSON results file from macsyfinder replaced by the TSV file for best solution

    my $gff = $repProkka."/".$RefSeq.".gff"; # the GFF file provided by prokka, example: NZ_AP017367/NZ_AP017367.gff
    if(-e $userGFF){
	$gff = $userGFF; # GFF file provided by user
    }
    else{
	$gff = $repProkka."/".$RefSeq.".gff";
    }

    my $resultsCRISPRs = $ResultDir."/Crisprs_REPORT.tsv";
    my $resultsTmpCRISPRs = $ResultDir."/TmpCrisprs_REPORT.tsv";

    open (RESULTS, ">>$results") or die "open : $!";  # Open $results in writing-adding mode

    if($writeFullReport){ # Modifying '-ccvr' option output - December 2017
    	open (CAS, ">>$allCas") or die "open : $!";
    }     

    if($logOption){
      print LOG "\n Writing in file Cas_REPORT.tsv\n"; #[$hour:$min:$sec] 
        #print LOG "[$hour:$min:$sec] Writing in file CRISPR-Cas_systems_vicinity.tsv\n";
    }

    print RESULTS "############################################\n";
    print RESULTS "$RefSeq ( $seqDesc ) \n";
    print RESULTS "--------------------------------------------\n";
    
    if( (-e $userGFF) and (-e $userFAA) ){}
    else{
	    if($logOption){
		if($useProdigal){
			print LOG "$prodigal\n"; #[$prokka_hour:$prokka_min:$prokka_sec] 
			#print "PRODIGAL command = $prodigal\n";
		}
		else{
	    		print LOG " $prokka\n"; #[$prokka_hour:$prokka_min:$prokka_sec]
		}
	    }
	    if($useProdigal){
	    	system($prodigal); # Launch Prodigal 
	    }
	    else{
		system($prokka); # Launch Prokka 
	    }
    }

    #Addition of a test to see if any CDS exist (DC, 2019 + modifications on may 2020)
    my $cdsCheck = 1;
    if (-z $proteome){ #Checking if $proteome file is empty
	if($logOption){
		print LOG " CDS file ($proteome) is empty! No Cas search because no CDS was detected!\n"; #[$macsyfinder_hour:$macsyfinder_min:$macsyfinder_sec]
	}
	if($quiet){}
	else{
	  	print "No Cas search because no CDS was detected!\n";
	}

	$cdsCheck = 0; #CasFinder will not be launched if no CDS have been detected ; $launchCasFinder has been replaced by $cdsCheck
	#next;
    }
    else{
	$cdsCheck = 1;
    }
    #End of Addition

    if (-e $gff and $cdsCheck){   # If GFF is created (Prokka is done) + addition of $cdsCheck
	if($logOption){
		print LOG " $macsyfinder\n"; #[$macsyfinder_hour:$macsyfinder_min:$macsyfinder_sec]
	}
	system($macsyfinder); # We can launch MacSyFinder
    }

    ## get the GFF file generated by Prodigal or Prokka (we will call it: 'annotation_$RefSeq.gff') 
    my $tmpGFFprokka = $ResultDir."/annotation_$RefSeq.gff"; #$ResultDir."/annotation_prokka_$RefSeq.gff";
    if(-e $gff and -d $ResultDir."/GFF"){

    	system ("cp $gff $tmpGFFprokka"); 

	makesystemcall("cp -r $tmpGFFprokka $ResultDir/GFF");
    }

    ## get the summary report which has been created by casfinder
    my $summaryCas = $casDir."best_solution_summary.tsv"; # the summary results file from macsyfinder
    if(-e $summaryCas and $gscf){
	# copy the file to Casfinder_summary.tsv
	my $tmpSummary = $ResultDir."/Casfinder_summary_$RefSeq.tsv";
	my $cmdCopySummary = "cp $summaryCas $tmpSummary";
	system($cmdCopySummary);
    }

    ### If JSON (TSV best solution result file) has been created    if ( $line =~ /$v/ && $line !~ /#/ ) {}

    if (-e $tsv){ 
		#my $p = json_file_to_perl ($json); #Translate JSON file into Perl object 

		#my $sizeTable = @{$p}; #Global size of main table of JSON file
		my $sizeTable = 0;
		
		my @tabGenes=(); # table containing sequences
		my %hSystem=();
		my %hGeneInfo=();
		my %hGeneSystem=();

		## Push all needed genes in a table
		# for (my $i = 0; $i < $sizeTable; $i++) {
			# my $sizeGenes2 = @{$p->[$i]->{'genes'}};
	
			# for (my $j = 0; $j < $sizeGenes2; $j++) {
				# my $sequenceID = $p->[$i]->{'genes'}->[$j]->{'id'};
				# push (@tabGenes, $sequenceID);
			# }
		# }
		open (TSVMAC,"<$tsv") or die $!;
		while (my $lineTSV = <TSVMAC>) {
			chomp($lineTSV);
			#next if($lineTSV =~ s/^\s*$//);
			if ($lineTSV !~ /#/ && $lineTSV !~ /^replicon\t/ && $lineTSV !~ /^\s*$/){	
			    #next if /^\s*$/;
			    my @tabTSV = split(/\t/, $lineTSV);	
				my $sequenceID = $tabTSV[1]; #index corresponding to gene ID
				push (@tabGenes, $sequenceID);
				my @tabTMPSys = split('_', $tabTSV[5]); # temp table to extract the system number
				$hSystem{$tabTMPSys[$#tabTMPSys]} = $tabTSV[5];
				#print "HSYSTEM value for System number: $tabTMPSys[$#tabTMPSys] \n";
				push (@{ $hGeneSystem{$tabTMPSys[$#tabTMPSys]} }, $tabTSV[1]);  # replacing $tabTSV[6] by $tabTMPSys[$#tabTMPSys]
				$hGeneInfo{$tabTSV[1]} = ["$tabTSV[2]", "$tabTSV[4]", "$tabTSV[6]", "$tabTSV[12]"]; #hash contains: gene name, systemCAS, sys_loci, hit_status
			}
		}
		$sizeTable = keys %hSystem;
		#print "SIZE TABLE IS: $sizeTable \n";
		close(TSVMAC);

		open (GFF, "<$gff") or die "open : $!";  #Open $gff in reading mode

		my $line=""; # To read file

		my %hashGeneType=(); #Hash table for sequence Type
		my %hashGeneBegin=(); #Hash table for sequence Begin
		my %hashGeneEnd=(); #Hash table for sequence End
		my %hashGeneStrand=(); #Hash table for sequence Strand
		my %hashGeneOther=(); #Hash table for sequence Other

			while (defined ($line = <GFF>)){
				foreach my $v (@tabGenes) {
				chomp($line);
				#print "MY GENE IS: $v , \n";
				my @partIDs = split(/_/, $v);
				my $elemIDpart = "ID=1_".$partIDs[$#partIDs].";";
				    if($useProdigal){
					if ( $line =~ /$elemIDpart/ && $line !~ /#/ ) {
			
						my @tab = split(/\t/, $line); #table containing wanted values of GFF file

						$hashGeneType{$v} = $tab[2]; #Type Seq
						$hashGeneBegin{$v} = $tab[3]; #Begin Seq
						$hashGeneEnd{$v} = $tab[4]; #End Seq
						$hashGeneStrand{$v} = $tab[6]; #Strand (orientation) of Seq
						$hashGeneOther{$v} = $tab[8]; #Other information of Seq

					}
				    }
				    else{
					# using prokka
					if ( $line =~ /$v/ && $line !~ /#/ ) {
			
						my @tab = split(/\t/, $line); #table containing wanted values of GFF file

						$hashGeneType{$v} = $tab[2]; #Type Seq
						$hashGeneBegin{$v} = $tab[3]; #Begin Seq
						$hashGeneEnd{$v} = $tab[4]; #End Seq
						$hashGeneStrand{$v} = $tab[6]; #Strand (orientation) of Seq
						$hashGeneOther{$v} = $tab[8]; #Other information of Seq

					}
				    }					
				}
			}

		### JSON Cas file For Cas systems
		my $jsonLineCas = "";
  		
		#open(JSONCAS,">",$ResultDir."/".$RefSeq."_cas.json") or die("Could not open the JSON file $RefSeq _cas.json because $!\n");
  		#print JSONCAS "[";
 		$jsonCAS = "[";
		
		#my $countSys = 1; # to count Cas systems in JSON file
		my $countGeneralCas = 0;

		#foreach my $i (keys %hSystem) {
		for (my $i = 1; $i <= $sizeTable; $i++) {   # loop for each system
			#print "############################################\n";
			if(! $quiet){
            		#print "### System: $p->[$i]->{'name'} ($p->[$i]->{'id'})\n";
					print "### System: $hSystem{$i}\n";
            		print "#SequenceID\tCas-type/subtype\tGene status\tSystem\tType\tBegin\tEnd\tStrand\tOther_information\n";
			}

			print RESULTS "### System: $hSystem{$i}\n";
			print RESULTS "#SequenceID\tCas-type/subtype\tGene status\tSystem\tType\tBegin\tEnd\tStrand\tOther_information\n";
			
			my $systemName = $hSystem{$i};
			my $stringCasGenes = "[";

			my @tabPositions = ();

			my $sizeGenes = @{ $hGeneSystem{$i} };
			
			my $secondChainCas = "\"Genes\": [\n"; # Cas genes for JSON file

			for (my $j = 0; $j < $sizeGenes; $j++) {
		
			    my $seqName = $hGeneSystem{$i}[$j]; #my $seqName = $p->[$i]->{'genes'}->[$j]->{'id'};
			    my $matcher = $hGeneInfo{$seqName}[0]; #$p->[$i]->{'genes'}->[$j]->{'match'};
 			    if ( defined($matcher) ){
				if(! $quiet){	
				print "$seqName\t";#"$p->[$i]->{'genes'}->[$j]->{'id'}\t";  #Retrieve genes 'id'
				print "$matcher\t";#"$p->[$i]->{'genes'}->[$j]->{'match'}\t"; #Retrieve genes 'match' corresponding to Cas type
				print "$hGeneInfo{$seqName}[3]\t"; #"$p->[$i]->{'genes'}->[$j]->{'gene_status'}\t"; #Retrieve genes 'status' corresponding to information 'mandatory', 'forbidden', 'accessory'
				print "$hGeneInfo{$seqName}[1]\t"; #"$p->[$i]->{'genes'}->[$j]->{'system'}\t"; #Retrieve genes 'system' corresponding to information on Cas types 'Type II, Type III, ...'
				print "$hashGeneType{$seqName}\t"; #"$tab[2]\t";
				print "$hashGeneBegin{$seqName}\t"; #"$tab[3]\t";
				print "$hashGeneEnd{$seqName}\t"; #"$tab[4]\t";
				print "$hashGeneStrand{$seqName}\t"; #"$tab[6]\t";
				print "$hashGeneOther{$seqName}\n"; #"$tab[8]\n";
				}
				#Print Same thing in RESULTS file

				print RESULTS "$seqName\t";  #Retrieve genes 'id'
				print RESULTS "$matcher\t"; #Retrieve genes 'match' corresponding to Cas type
				print RESULTS "$hGeneInfo{$seqName}[3]\t"; #Retrieve genes 'status' : 'mandatory', 'forbidden', 'accessory'
				print RESULTS "$hGeneInfo{$seqName}[1]\t"; #Retrieve genes 'system' or Cas types: 'Type II, Type III, ...'
				print RESULTS "$hashGeneType{$seqName}\t";
				print RESULTS "$hashGeneBegin{$seqName}\t";
				print RESULTS "$hashGeneEnd{$seqName}\t";
				print RESULTS "$hashGeneStrand{$seqName}\t"; 
				print RESULTS "$hashGeneOther{$seqName}\n"; 

				#Create Summary System data 				
				$stringCasGenes .= "$matcher ($hashGeneBegin{$seqName},$hashGeneEnd{$seqName},$hashGeneStrand{$seqName}); ";

				#Fill Cas genes info into JSON Cas file
				$secondChainCas .= "{\n\"Sub_type\": \"".$matcher."\",\"Start\": ".$hashGeneBegin{$seqName}.",\"End\": ".$hashGeneEnd{$seqName}.",\"Orientation\": \"".$hashGeneStrand{$seqName}."\"\n},";

				#Push positions in table
				push(@tabPositions,$hashGeneBegin{$seqName});
				push(@tabPositions,$hashGeneEnd{$seqName});


			       ## Write raw Cas genes sequences
			       # get subsequence from genes Fasta file created by Prokka
			       
			       my $inputfileSeq = $repProkka."/".$RefSeq.".ffn"; # the .ffn file provided by prokka, example: NZ_AP017367/NZ_AP017367.ffn
			       if(-e $inputfileSeq){
					my $casSeqFile = $ResultDir."/rawCas.fna";
					open (RAWCAS, ">>$casSeqFile") or die "open : $!";  #Open casSeqFile

				      	my $dbSeq = Bio::DB::Fasta->new($inputfileSeq); # sequence file (or $userfile) 
				       
				       	# DC - 11/2017 - get raw sequence of Cas
				       	my $rawCASsequence = $dbSeq->seq($seqName); # Cas gene sequence
					
					my $len = 80;

					my $formatted_seq = ">".$RefSeq."|".$seqName."|".$systemName."|".$matcher." ".$hashGeneBegin{$seqName}.",".$hashGeneEnd{$seqName}."\n";
					while (my $chunk = substr($rawCASsequence, 0, $len, "")) {
						$formatted_seq .= "$chunk\n";
					}

					print RAWCAS $formatted_seq; # ."\n";

					close(RAWCAS);
			       }
			       else{
					$inputfileSeq = $RefSeq.".fna";
				     if(-e $inputfileSeq){
					my $casSeqFile = $ResultDir."/rawCas.fna";
					open (RAWCAS, ">>$casSeqFile") or die "open : $!";  #Open casSeqFile

				      	my $dbSeq = Bio::DB::Fasta->new($inputfileSeq); # sequence file (or $userfile) 
				        my $seqidCas = $seq->id; # ID of sequence
				       	# DC - 11/2017 - get raw sequence of Cas
				       	my $rawCASsequence = $dbSeq->seq($seqidCas, $hashGeneBegin{$seqName} => $hashGeneEnd{$seqName});#$dbSeq->seq($seqName); # Cas gene sequence
					####
					#my $rawCRISPRsequence = $dbSeq->seq($seqid, $start => $end);
					####
					my $len = 80;

					my $formatted_seq = ">".$RefSeq."|".$seqName."|".$systemName."|".$matcher." ".$hashGeneBegin{$seqName}.",".$hashGeneEnd{$seqName}."\n";
					while (my $chunk = substr($rawCASsequence, 0, $len, "")) {
						$formatted_seq .= "$chunk\n";
					}

					print RAWCAS $formatted_seq; # ."\n";

					close(RAWCAS);
   				     }  
			       }
			      ##

			    }

			} ## End FOR genes loop
			my $beginCasCluster = min(@tabPositions);
			my $endCasCluster = max(@tabPositions);

			$stringCasGenes .= "]";
			$stringCasGenes =~ s/; ]/]/;

			if(! $quiet){print "####Summary system $hSystem{$i}:begin=$beginCasCluster;end=$endCasCluster;sequenceID=$RefSeq\n\n";}
			print RESULTS "####Summary system $hSystem{$i}:begin=$beginCasCluster;end=$endCasCluster:{sequenceID=$RefSeq} : $stringCasGenes\n\n";

			## Fill CCS file and count Cas systems
			my $ccsCas = $systemName."[".$beginCasCluster.";".$endCasCluster."],";
            push (@tabCCScas, $ccsCas);
            $hashCountCas{$systemName}++;
            ##

			### Fill Cas Table with name of system (if $clusteringThreshold is set)
			if($clusteringThreshold){
				
				#if($p->[$i]->{'name'} eq "General-CAS"){
				$countGeneralCas++;
				my $otherNameGeneralCas = $hSystem{$i}."_n".$countGeneralCas;
				$hashCasBegin{$otherNameGeneralCas} = $beginCasCluster;
	  			$hashCasEnd{$otherNameGeneralCas} = $endCasCluster;	
				#}
				#$hashCasBegin{$p->[$i]->{'name'}} = $beginCasCluster;
	  			#$hashCasEnd{$p->[$i]->{'name'}} = $endCasCluster;
				
				push (@otherTabCas, $otherNameGeneralCas); #General-CAS
			}
			###


			## Fill the JSON file dedicated to Cas systems
			$secondChainCas .= "]},";
			$secondChainCas =~ s/,]},/]},/;
			$jsonLineCas .= "{\n\"Type\": \"".$hSystem{$i}."\",\"Start\": ".$beginCasCluster.",\"End\": ".$endCasCluster.",".$secondChainCas; 
			## End Fill JSON Cas

			# Read file containing CRISPRs and watch if one is close to a Cas systems
			# Get also the Consensus repeat (or DR)
			@otherTabCrispr = ();
			open (CRISPR, "<$resultsCRISPRs") or die "open : $!";
			#open (TMP, ">>$resultsTmpCRISPRs") or die "open : $!"; #$resultsTmpCRISPRs
			while(<CRISPR>) {
    				chomp($_);
				#if ($_ =~  m/Strain\tSequence/) {
				#    print TMP "$_\n";
				#}
				if ($_ =~  m/$RefSeq/) {
        
				    my @tab = split(/\t/, $_);
				    if(defined($tab[5])){	# December 2017 - addition of and $writeFullReport option

					### Group all CRISPRs in a table (if $clusteringThreshold is set)
					if($clusteringThreshold){
						push (@otherTabCrispr, $tab[4]);
						$hashCrisprBegin{$tab[4]} = $tab[5];
	  					$hashCrisprEnd{$tab[4]} = $tab[6];
					}
					###
					if( (($tab[5]-$vicinity) <= $endCasCluster) and (($tab[5]-$vicinity) > $beginCasCluster) and $writeFullReport ){
						#$_ .= "\t$p->[$i]->{'name'} [$beginCasCluster;$endCasCluster]\n";
						#print TMP $_;
						#print CRISPR "\t$p->[$i]->{'name'} [$beginCasCluster;$endCasCluster]\n";
						print CAS "$tab[4]\t$tab[4] [$tab[5];$tab[6]]\t";
						print CAS "$hSystem{$i} [$beginCasCluster;$endCasCluster]\n";
					}
					elsif ( (($tab[6]+$vicinity) >= $beginCasCluster) and (($tab[6]+$vicinity) < $endCasCluster) and $writeFullReport ){
						#$_ .= "\t$p->[$i]->{'name'} [$beginCasCluster;$endCasCluster]\n";
						#print TMP $_;
						#print CRISPR "\t$p->[$i]->{'name'} [$beginCasCluster;$endCasCluster]\n";
						print CAS "$tab[4]\t$tab[4] [$tab[5];$tab[6]]\t";
						print CAS "$hSystem{$i} [$beginCasCluster;$endCasCluster]\n";
					} 
				
				    }
				}
			}
			#print CAS "$p->[$i]->{'name'}\t$p->[$i]->{'id'}\t$RefSeq\t$beginCasCluster\t$endCasCluster\n";
		} # end of for loop on table of Cas system
		print "\n";
		print RESULTS "\n";
		
		$nbCas = @tabCCScas;
		my $chainCountCas = "";
		if(%hashCountCas){
		    foreach my $elemCountCas (sort keys %hashCountCas) {
		        $chainCountCas .= $elemCountCas." (n=".$hashCountCas{$elemCountCas}."), ";
		    }
		}
		print RESULTCCS "@tabCCScas\t$nbCas\t$chainCountCas\n";
		close (RESULTCCS); # Close RESULTS file

		close (GFF); # Close GFF file 
		close (RESULTS); # Close RESULTS file
		if($writeFullReport){
			close (CAS); # Close CAS file
		}
		close (CRISPR);
		#close (TMP);
		#makesystemcall("mv $resultsTmpCRISPRs $resultsCRISPRs"); # Move TMP file into CRISPR

		##JSON Cas (note replacement of '}]' by '}')
		$jsonLineCas .= "]\n";
  		$jsonLineCas =~ s/,]/]/; # change end of file to better fit JSON format
  
      		$jsonLineCas =~ s/,/,\n/g; # make file more readable

  		#print JSONCAS $jsonLineCas;
		$jsonCAS .= $jsonLineCas;

  		#close JSONCAS;
		## End JSON Cas

		### Generate clusters of CRISPR-Cas in function of threshold
		if($clusteringThreshold){
			my @otherTabCCC = (); # table to store CRISPR and Cas in the good order
			my $countClust = 0;
			my $nTabCrispr = @otherTabCrispr; # Nb of elements in @otherTabCrispr
			my $nTabCas = @otherTabCas; # same thing for @otherTabCas

			my ($i, $j, $k) = (0, 0, 0);

			while($i < $nTabCrispr and $j < $nTabCas){
				if( $hashCrisprBegin{$otherTabCrispr[$i]} < $hashCasBegin{$otherTabCas[$j]} ){
				  push (@otherTabCCC, $otherTabCrispr[$i]);
				  $i++;
				}
				else{
				  push (@otherTabCCC, $otherTabCas[$j]);
				  $j++;
				}
			}

			while($i < $nTabCrispr){
				push (@otherTabCCC, $otherTabCrispr[$i]);
				$i++;
			}

			while($j < $nTabCas){
				push (@otherTabCCC, $otherTabCas[$j]);
				$j++;
			}
			#my $nbInTabCCC = @otherTabCCC;
			#print "Nb CRISPR = $nTabCrispr ; Nb Cas = $nTabCas ;\n TAB CCC (nb:$nbInTabCCC) = @otherTabCCC\n\n";
			
			#Fill summary table with essantial info
			foreach my $elemVccc (@otherTabCCC) {
				my $tmpStatus = "";
				if( defined($hashCrisprBegin{$elemVccc})  ){

					$tmpStatus = $elemVccc."[".$hashCrisprBegin{$elemVccc}.";"
							.$hashCrisprEnd{$elemVccc}."]"; 

					push(@tabCRISPRCasClusters, $tmpStatus);
				}
				elsif( defined($hashCasBegin{$elemVccc}) ){
					$tmpStatus = $elemVccc."[".$hashCasBegin{$elemVccc}.";"
							.$hashCasEnd{$elemVccc}."]"; 
					push(@tabCRISPRCasClusters, $tmpStatus);
				}
			}
			#print "SUMMARY TABLE= @tabCRISPRCasClusters \n";
			### get clusters of CRISPRs or Cas in function of threshold
			my $tempCCC = "";
			my @tabCCC3 = ();
			my $varL = 0;
			my $exitValue = 0;
			my $clusteringCC = 1;
			$k=0;
			#if( defined($hashCrisprBegin{$otherTabCCC[$k]})  ){
			#	  $tempCCC = $otherTabCCC[$k]."[".$hashCrisprBegin{$otherTabCCC[$k]}.";"
			#				.$hashCrisprEnd{$otherTabCCC[$k]}."], "; 

			#}
			#elsif( defined($hashCasBegin{$otherTabCCC[$k]})){  # Same For Cas
			#	  $tempCCC = $otherTabCCC[$k]."[".$hashCasBegin{$otherTabCCC[$k]}.";"
			#				.$hashCasEnd{$otherTabCCC[$k]}."], "; 
				  ##
			#}
			#$k+=1;
            while ( ($k <= $#otherTabCCC) and ($varL <= $#otherTabCCC) ) {
                #print "11111111111\n";
				if( defined($hashCrisprBegin{$otherTabCCC[$k]})  ){
					  $tempCCC = $otherTabCCC[$k]."[".$hashCrisprBegin{$otherTabCCC[$k]}.";"
								.$hashCrisprEnd{$otherTabCCC[$k]}."],"; 

				}
				elsif( defined($hashCasBegin{$otherTabCCC[$k]})){  # Same For Cas
					  $tempCCC = $otherTabCCC[$k]."[".$hashCasBegin{$otherTabCCC[$k]}.";"
								.$hashCasEnd{$otherTabCCC[$k]}."],"; 
				
				}				
				$varL=$k+1;
				$clusteringCC = 1;
				while ( ($varL<=$#otherTabCCC) and $clusteringCC ) { #

				  
				  if( defined($hashCrisprBegin{$otherTabCCC[$k]}) and defined($hashCrisprBegin{$otherTabCCC[$varL]}) and
				      ($hashCrisprBegin{$otherTabCCC[$varL]} - $hashCrisprEnd{$otherTabCCC[$k]}) 
					< $clusteringThreshold) {

				    $tempCCC .= $otherTabCCC[$varL]."[".$hashCrisprBegin{$otherTabCCC[$varL]}.";"
						.$hashCrisprEnd{$otherTabCCC[$varL]}."],";

				    
                    #$exitValue++;
				    
				  }
			 	  elsif( defined($hashCasBegin{$otherTabCCC[$varL]}) and defined($hashCrisprEnd{$otherTabCCC[$k]}) and
					 ($hashCasBegin{$otherTabCCC[$varL]} - $hashCrisprEnd{$otherTabCCC[$k]}) 
					< $clusteringThreshold ){

				    $tempCCC .= $otherTabCCC[$varL]."[".$hashCasBegin{$otherTabCCC[$varL]}.";"
						.$hashCasEnd{$otherTabCCC[$varL]}."],";
				    #$exitValue++;
				  }
				  elsif( defined($hashCasBegin{$otherTabCCC[$k]}) and defined($hashCasBegin{$otherTabCCC[$varL]}) and
					 ($hashCasBegin{$otherTabCCC[$varL]} - $hashCasEnd{$otherTabCCC[$k]}) 
					< $clusteringThreshold ){

				    $tempCCC .= $otherTabCCC[$varL]."[".$hashCasBegin{$otherTabCCC[$varL]}.";"
						.$hashCasEnd{$otherTabCCC[$varL]}."],";
				    #$exitValue++;
				  }
				  elsif( defined($hashCasBegin{$otherTabCCC[$k]}) and defined($hashCrisprBegin{$otherTabCCC[$varL]}) and
					 ($hashCrisprBegin{$otherTabCCC[$varL]} - $hashCasEnd{$otherTabCCC[$k]}) 
					< $clusteringThreshold ){

				    $tempCCC .= $otherTabCCC[$varL]."[".$hashCrisprBegin{$otherTabCCC[$varL]}.";"
						.$hashCrisprEnd{$otherTabCCC[$varL]}."],";
				    #$exitValue++;
				  }


				  else{
					$clusteringCC = 0;
				  }

				  $k = $varL;
				  $varL++;
				  #$k+=1;
			  	}

			  chop($tempCCC);
			  push (@tabCCC3, $tempCCC);
			  #print "tempCCC (second loop) = $tempCCC ;\nvarL = $varL\nk = $k\nexitValue = $exitValue \n";
			  $tempCCC = "";

			}

						
			open (RESULTCCC, ">>$clusterResults") or die "open : $!";
			foreach my $vCCC (@tabCCC3) {
				my @tabTemp2 = split(/,/, $vCCC);
				#print "@tabTemp2\n";
				my ($startCCC, $endCCC, $nbCrCCC, $nbCasCCC);
                $nbCrCCC = 0;
                $nbCasCCC = 0;
				my $nbTmp = @tabTemp2;
				
				if($nbTmp > 1){
				my ($idElemCCfirst, $otherCC1) = split(/\[/,$tabTemp2[0]);
				my ($idElemCClast, $otherCC2) = split(/\[/,$tabTemp2[$#tabTemp2]);
				

				if( defined($hashCasBegin{$idElemCCfirst}) ){
					$startCCC = $hashCasBegin{$idElemCCfirst};
				}
				elsif( defined($hashCrisprBegin{$idElemCCfirst}) ){
					$startCCC = $hashCrisprBegin{$idElemCCfirst};
				}
                else{
                    $startCCC = 0;
                }
				
				if( defined($hashCasEnd{$idElemCClast}) ){
					$endCCC = $hashCasEnd{$idElemCClast};
				}
				elsif( defined($hashCrisprEnd{$idElemCClast}) ){
					$endCCC = $hashCrisprEnd{$idElemCClast};
				}
				
				foreach my $dataCCC (@tabTemp2) {
					my ($elemCCC, $otherCC2) = split(/\[/,$dataCCC);
					#print "$elemCCC\n";
					if( defined($hashCasBegin{$elemCCC}) ){
						$nbCasCCC++;
					}
					elsif( defined($hashCrisprBegin{$elemCCC}) ){
						$nbCrCCC++;
					}
				}
				
				$countClust++; 
				print "Cluster nb $countClust = $vCCC\n";
				print RESULTCCC "$RefSeq\t$countClust\t$nbCrCCC\t$nbCasCCC\t$startCCC\t$endCCC\t$vCCC\n";
				}
				
			}
			close(RESULTCCC);
			
		}
		### End generation of clusters

	} # End of IF (-e $JSONfile) { ... }
	else {
		if (! $quiet) { print "No Cas results\n"; }
		$nbCas = 0;
		print RESULTCCS "@tabCCScas\t0\n";
		close (RESULTCCS); # Close RESULTS file
	}

    ###

  };

  if($@){
    #An error occurred...
    print "An error occurred in CasFinder function\n";
  }
  #print "Final NB CAS = $nbCas\n";

  #create directories Prokka and CasFinder then move all related files to this directory (if -cas is set)
  ##if ($launchCasFinder){
	if($userGFF and $userFAA){}
	else{
	    if($useProkka){
		my $directoryProkka = $ResultDir."/Prokka";  # Directory for Prokka
		mkdir $directoryProkka unless -d $directoryProkka;
		if(-d $repProkka and -e $repProkka and -e $directoryProkka){
			makesystemcall("cp -r $repProkka $directoryProkka");
	        }
	    }
	    else{
		my $directoryProdigal = $ResultDir."/Prodigal";  # Directory for Prodigal
		mkdir $directoryProdigal unless -d $directoryProdigal;

		if(-d $repProkka){
			makesystemcall("cp -r $repProkka $directoryProdigal");
		}
	    }
	}

	my $directoryCasfinder = $ResultDir."/CasFinder";  # Directory for CasFinder
	mkdir $directoryCasfinder unless -d $directoryCasfinder;
	if(-d $casDir){
	makesystemcall("cp -r $casDir $directoryCasfinder");
	}
  ##}

  return ($nbCas, $results, $jsonCAS, @tabCRISPRCasClusters);
  
}
#------------------------------------------------------------------------------
# DC - 06/2017 - Rewrite Report adding Cas genes information
sub fullReport{
  my($ResultDir)=@_;
  my $resultsCRISPRs = $ResultDir."/Crisprs_REPORT.tsv";
  my $allCas = $ResultDir."/CRISPR-Cas_systems_vicinity.tsv";
  my $fullReport = $ResultDir."/CRISPR-CAS_vicinity_REPORT.tsv";
  my %hashCrisprCas = ();

  if( (-e $allCas) and (-e $resultsCRISPRs) ){
	  # Read file containing all Cas and save data
	  open (CAS, "<$allCas") or die "open : $!";
	  while(<CAS>) {
	  	chomp($_);
		my ($crisprID,$crispr_positions,$cas) = split(/\t/, $_);
		#print "ID= $crisprID\n";
	 	if(defined($crisprID)){
			$hashCrisprCas{$crisprID} = $cas;
			#print "CAS= $hashCrisprCas{$crisprID}\n";
		}
	  }
	  close (CAS);

	  # Read file containing CRISPRs and write fullReport
	  open (CRISPR, "<$resultsCRISPRs") or die "open : $!";
	  open (FULL, ">$fullReport") or die "open : $!";
	  #my ($hour,$min,$sec) = Now();

	  if($logOption){
	  	print LOG "\n Create CRISPR-CAS_vicinity_REPORT.tsv (corresponding to CRISPRs arrays and neighboring CAS systems in a range of $vicinity bps)\n";
	  } #[$hour:$min:$sec]

	  while(<CRISPR>) {
	  	chomp($_);
		
		my @tab = split(/\t/, $_);

		if($_ =~  m/Strain\tSequence/){
			print FULL "$_\tNeighbor_Cas (vicinity=$vicinity bps)\n";
		}
		elsif ( (defined($tab[4])) and ($_ =~  m/$tab[4]/) and (defined($hashCrisprCas{$tab[4]})) ) {
			print FULL "$_\t$hashCrisprCas{$tab[4]}\n";
			#print $hashCrisprCas{$tab[4]}."\n";
		}
		else{
			print FULL "$_"."\n";
		}
	  }

	  close (CRISPR);
	  close (FULL);
  }
  makesystemcall("cp $fullReport $resultsCRISPRs"); # Move fullReport file into CRISPRs_report
  #makesystemcall("rm -f $allCas"); # remove $allCas
  return $fullReport;
}
#------------------------------------------------------------------------------
#DC - 06/2017 - Count CRISPRs orientation
sub countOrientation{
  my($ResultDir,$nbcrispr)=@_;
  my $resultsCRISPRs = $ResultDir."/Crisprs_REPORT.tsv";
  my $orientationCountsFile = $ResultDir."/crisprs_orientations_count.tsv";

  my $countMatchingOrientation = 0;
  my $countForwardFinder = 0;
  my $countReverseFinder = 0;
  my $countForwardCRISPRDirection = 0;
  my $countReverseCRISPRDirection = 0;


  if( -e $resultsCRISPRs ){

	# Read file containing CRISPRs and write CRISPRs orientations counts (based on repeats)
	open (CRISPR, "<$resultsCRISPRs") or die "open : $!";
		
	while(<CRISPR>) {
	  	chomp($_);
		
	    my @tab = split(/\t/, $_);
	    my $potentialDir = $tab[8];
	    my $crisprDir = $tab[9];
	    if ( (defined($potentialDir)) or (defined($crisprDir)) ){	

		if($potentialDir =~  m/^F/){
			$countForwardFinder++;
		}
		elsif($potentialDir =~  m/^R/){
			$countReverseFinder++;
		}

		if($crisprDir =~  m/^F/){
			$countForwardCRISPRDirection++;
		}
		elsif($crisprDir =~  m/^R/){
			$countReverseCRISPRDirection++;
		}
		
		if( ($potentialDir =~  m/^F/) and ($crisprDir =~  m/^F/) ){
			$countMatchingOrientation++;
		}
		elsif(($potentialDir =~  m/^R/) and ($crisprDir =~  m/^R/)){
			$countMatchingOrientation++;
		}
	     }
	}

	close (CRISPR);
  }

  open (ORI, ">$orientationCountsFile") or die "open : $!";
  #my ($hour,$min,$sec) = Now();

  if($logOption){
  	print LOG "\n Create orientationCounts file ($orientationCountsFile)\n"; #[$hour:$min:$sec]
  }

  print ORI "#### Statistics on CRISPRs orientation by CRISPRCasFinder vs. CRISPRDirection\n\n";
  print ORI "Total number of CRISPRs arrays found = $nbcrispr\n";
  print ORI "Number of perfect macthes between CRISPRCasFinder's potential orientation (based on AT%) and CRISPRDirection = $countMatchingOrientation\n";
  print ORI "Number of Forward by CRISPRCasFinder = $countForwardFinder\n";
  print ORI "Number of Forward by CRISPRDirection = $countForwardCRISPRDirection\n";
  print ORI "Number of Reverse by CRISPRCasFinder = $countReverseFinder\n";
  print ORI "Number of Reverse by CRISPRDirection = $countReverseCRISPRDirection\n\n";

  print ORI "Number of unoriented by CRISPRCasFinder = ".($nbcrispr - ($countForwardFinder + $countReverseFinder))."\n";
  print ORI "Number of unoriented by CRISPRDirection = ".($nbcrispr - ($countForwardCRISPRDirection + $countReverseCRISPRDirection))."\n";

  if(! $quiet){
  print "#### Statistics on CRISPRs orientation by CRISPRCasFinder vs. CRISPRDirection\n\n";
  print "Total number of CRISPRs arrays found = $nbcrispr\n";
  print "Number of perfect macthes between CRISPRCasFinder and CRISPRDirection = $countMatchingOrientation\n";
  print "Number of Forward by CRISPRCasFinder = $countForwardFinder\n";
  print "Number of Forward by CRISPRDirection = $countForwardCRISPRDirection\n";
  print "Number of Reverse by CRISPRCasFinder = $countReverseFinder\n";
  print "Number of Reverse by CRISPRDirection = $countReverseCRISPRDirection\n\n";

  print "Number of unoriented by CRISPRCasFinder = ".($nbcrispr - ($countForwardFinder + $countReverseFinder))."\n";
  print "Number of unoriented by CRISPRDirection = ".($nbcrispr - ($countForwardCRISPRDirection + $countReverseCRISPRDirection))."\n";
  }
  return $orientationCountsFile;
}
#------------------------------------------------------------------------------
#LK
# parse result files and create a GFF file for all
sub makeGff{
  my($ResultDir,$inputfile,$seqDesc,$nbrcris,$OneSpacerCris_nbr)=@_;
  my @idDir = (); # DC - table containing crisprs IDs and potential directions
  my @fastaExtensions=qw(.fasta .fna .faa .mfa .fa .fas .ffn);
  my $inputBasename=basename($inputfile,@fastaExtensions);
  my @dir=<$ResultDir/$inputBasename*>;
  my $globalAT2 = "";

  #my ($hour,$min,$sec) = Now();

  if($logOption){
  	print LOG "\n Create GFF3 file ($inputBasename.gff) corresponding to CRISPRs arrays\n"; #[$hour:$min:$sec] 
  }
  # parse
  my $GFFstr="";
  open GFF,">$ResultDir/$inputBasename.gff" or die "Could not open output GFF in $ResultDir/ because $!\nI will not create the GFF file.\n";
  for my $dir(@dir){
    my @spacerFile=<$dir/Spacers_*>;
    for my $spacerFile(@spacerFile){
      my($spacerNumber);
      if($spacerFile=~/_(\d+)$/){
        $spacerNumber=$1 ;
      } else {
        warn "Warning: Cannot understand result file $spacerFile. Skipping.\n";
        next;
      }

      # transform the report file to a GFF line
      my $reportFile="$dir/$inputBasename"."_PossibleCrispr_$spacerNumber";
      $reportFile="$dir/$inputBasename"."_Crispr_$spacerNumber" unless -e $reportFile;
      unless(-e $reportFile){
         warn "Impossible to find details for spacer $spacerNumber";
	 next;
      }
      my ($feature,$crisprID,$potentialDir,$globalAT)=reportToGff($reportFile,$inputfile);

      $globalAT2 = $globalAT; # to use globalAT value outside

      my $id_Dir = $crisprID."_".$potentialDir; # DC - to retrieve potential direction of CRISPRs
      push (@idDir,$id_Dir); # DC

      # add onto the GFF
      $GFFstr.="$feature\n" if defined $feature;
    }
  }

  $GFFstr="##gff-version 3; $nbrcris CRISPR(s); strain: $seqDesc \n$GFFstr"; # DC
  print GFF $GFFstr;
  close GFF;
  #print STDERR "GFF results are located in $ResultDir/$inputBasename.gff\n";

  return ("$ResultDir/$inputBasename.gff",@idDir);
}
#------------------------------------------------------------------------------
#LK
# transform a given report file to a line for a GFF file
# returns a string TODO a feature object
# http://www.bioperl.org/wiki/GFF3
# http://doc.bioperl.org/releases/bioperl-current/bioperl-live/Bio/SeqFeature/Generic.html
sub reportToGff{
  my($reportFile,$inputfile)=@_;
  #GFF columns
  my($seqid,$start,$end,$strand,$attributes,%tag);
  my($source,$type,$score,$phase)=("CRISPRCasFinder","CRISPR",".",".");
  $type="CRISPR" if $reportFile =~ /PossibleCrispr_\d+$/;  # DC replaced $type="PossibleCRISPR" by $type="CRISPR"
  
  my($featureId,$gffLine,@spacerInfo,$DRsequence);
  return unless -e $reportFile;
  # get data from the report file
  my $onSpacersList = 0;
  my $potentiallyFalse = 0; # DC - indicate if a given CRISPR is potentially false
  my $potentiallyFalseCRISPR = ""; # DC - corresponding text
  open(REPORT,"<",$reportFile) or die("Could not open the report file $reportFile because $!\n");

  while(<REPORT>){
    my($spacerStart,$spacerEnd,$spacerSequence);

    if(/#Potentially false/){ # DC - control if the given CRISPR is potentially false (#Potentially )
      $potentiallyFalse = 1; # boolean value to indicate false CRISPR (in case of duplicated spacers or too long/short spacers)
    }

    if(/Id: (\S+)/){
      $seqid=$1;

    }
    elsif (/Crispr.+?(\d+)\s+Crispr\D+(\d+)/){
      $start=$1;
      $end=$2;
      my $nextLine=<REPORT>;
      $nextLine=~s/://g; # remove colons
      my @field=split(/\s+/,$nextLine);
      shift(@field); # remove the hash sign
      for(my $i=0;$i<@field;$i+=2){
        $attributes.="$field[$i]=$field[$i+1];";
        $tag{$field[$i]}=$field[$i+1];
      }
    }
    elsif (/Spacer_begin_position/){
      $onSpacersList = 1;
    }
    elsif ($onSpacersList){
      next if /^#/;
      next if /^Spacers divisions:/;  # DC, to avoid info on long spacers
      next if /^\s*$/;
      my @field=split(/\s+/,$_);
      shift(@field);
      my $name=join("_","spacer",$field[0],$field[1]);
      push(@spacerInfo,{
        start=>$field[0],end=>($field[0]+$field[1]-1),
        type=>"CRISPRspacer",
        attributes=>"sequence=$field[2];Name=$name",
        }
      );
    }

  }
  close REPORT;
  
  # DC - 05/2017 - indicate false CRISPR
  if($potentiallyFalse == 1){
    $potentiallyFalseCRISPR = "potentially_false_CRISPR";
  } else{
    $potentiallyFalseCRISPR = "";
  }

  # derivations of the data - removed by DC - 05/2017
  #if($start<$end){
  #  $strand='1';
  #} else{
  #  $strand='0';
  #}

  # DC - 05/2017 - replace "//" and get real ID 
  my $reportDR = $reportFile; #DC
  $reportDR =~ s/\/\//\//; #DC - replace '//' by '/'
  my @tabReport = split(/\//, $reportDR); # DC
  my $idCRISPR = pop(@tabReport); # DC - get last element of table
  $featureId=join("_",$seqid,$start,$end);
  #$attributes.="name=$featureId;ID=$featureId";
  $attributes.="Name=$featureId;ID=$idCRISPR"; # DC - 05/2017
  $DRsequence=$tag{DR}; delete($tag{DR});
  
  #my ($hour,$min,$sec) = Now();

  if($logOption){
  	print LOG "\n --> GFF file ($inputfile.gff) with CRISPR ID = $idCRISPR and Consensus DR = $DRsequence\n"; #[$hour:$min:$sec]
  }


  my @tabIdCRISPR = split(/_/, $idCRISPR); # DC

  # DC adding DR sequence in GFF report (05/2017)
  my $dbSeq = Bio::DB::Fasta->new($inputfile); # DC- to store the fasta sequence (or $userfile)
  #my @idsSeq = $dbSeq->get_all_primary_ids; # DC- to get ID
  my $sequenceDR = ""; #DC - DR sequence
  my $globalSeq = $dbSeq->seq($seqid); # whole sequence to calculate global AT%

  # DC - 05/2017 - adding flanking sequence (left)
  my $leftFlankSeq = "";

  if( ($start-$flankingRegion) > 0 ){
  	$leftFlankSeq = $dbSeq->seq($seqid, ($start-$flankingRegion) => ($start-1)); # DC - Left flanking sequence
  }
  else{
	$leftFlankSeq = "UNKNOWN";
  }

  my $globalAT = atpercent($globalSeq); # Calculate AT% of the entire sequence

  # DC - 05/2017 - Calculate AT% in left flank
  my $leftAT = atpercent($leftFlankSeq); # DC - using function atpercent
  my $flankAttributes = "sequence=$leftFlankSeq;at%=".$leftAT; # DC

  # DC - 05/2017 - same thing for flanking sequence (right)
  my $rightFlankSeq = "";
  my $endGlobalSeq = length($globalSeq); # total length of sequence
  
  if( $end < ($endGlobalSeq-$flankingRegion) ){
  	$rightFlankSeq = $dbSeq->seq($seqid, ($end+1) => ($end+$flankingRegion)); # DC
  }
  else{
	$rightFlankSeq = "UNKNOWN";
  }

  my $rightAT = atpercent($rightFlankSeq); # DC - using function atpercent
  my $leader = ""; # DC - first leader corresponding to left flanking sequence
  my $leader2 = ""; # DC - second leader corresponding to right flanking sequence

  ## Repeat ID + Expert Orientation
  my @tabR = repeatIDandNb($DRsequence);
  my $expertOrientation = $tabR[2];
  
  my $potentialDir = ""; #DC - potential CRISPR orientation

  if( defined($expertOrientation) and ($expertOrientation =~  m/^F/) ){
	$potentialDir = "Forward";
      	$leader = ";leader";
      	$strand = "+";
  }
  elsif(defined($expertOrientation) and ($expertOrientation =~  m/^R/) ){
	$potentialDir = "Reverse";
      	$leader2 = ";leader";
      	$strand = "-";
  }
  elsif($leftAT == 0 || $rightAT == 0 || $leftAT == $rightAT){
    $potentialDir = "Unknown";
    $strand = ".";
  }
  else {
    if( ($leftAT > $rightAT) and ($leftAT > $globalAT) ){
      $potentialDir = "Forward";
      $leader = ";leader";
      $strand = "+";
    }
    elsif( ($leftAT < $rightAT) and ($rightAT > $globalAT) ) {
      $potentialDir = "Reverse";
      $leader2 = ";leader";
      $strand = "-";
    }
    else{
      $potentialDir = "Unknown";
      $strand = ".";
    }
  }


  # DC - Left Flanking Sequence
  $gffLine=join("\t",$seqid,$source,"LeftFLANK",($start-$flankingRegion),($start-1),$score,$strand,$phase,$flankAttributes.$leader.";Parent=$featureId;ID=fl_$start")."\n"; # DC - show results

  # Construct the gff line with 1) the entire CRISPR and then 2)DRs interspersed with spacers
  $gffLine.=join("\t",$seqid,$source,$type,$start,$end,$score,$strand,$phase,$attributes.";potential_direction=".$potentialDir.";".$potentiallyFalseCRISPR)."\n";
  my $DRstart=$start;
  
  # DC adding fasta file containing DRs 'DRLIST' (05/2017)
  #my $drFile = $reportDR."_DR.fna"; # DC
  #print "REPORT DR : $reportDR\n"; # DC
  my $idCrispr = pop(@tabIdCRISPR);
  my $drFile = $tabReport[0]."/".$tabReport[1]."/"."DRs_".$idCrispr; # DC
  

  open(DRLIST,">",$drFile) or die("Could not write file $drFile because $!\n"); # DC
  my $chaineDR = ""; # DC
 
  for my $s (@spacerInfo){

    $sequenceDR =  $dbSeq->seq($seqid, $DRstart => ($DRstart+$tag{DR_length}-1)); # DC - DR sequence

    # DRs surround spacers.

    # DC - 05/2017    
    #$gffLine.=join("\t",$seqid,$source,"CRISPRdr",$DRstart,($DRstart+$tag{DR_length}-1),$score,$strand,$phase,"Parent=$featureId")."\n";
    $gffLine.=join("\t",$seqid,$source,"CRISPRdr",$DRstart,($DRstart+$tag{DR_length}-1),$score,$strand,$phase,"sequence=$sequenceDR;Parent=$featureId;ID=DR_$DRstart")."\n"; # DC - 05/2017 Adding DR sequence
    
    # DC - 05/2017 - print DRs into fasta file
    $chaineDR = ">DR_".$DRstart."_".($DRstart+$tag{DR_length}-1); #DC
    print DRLIST "$chaineDR\n$sequenceDR\n"; #DC

    # spacer
    $$s{attributes}.=";Parent=$featureId;ID=sp_".$$s{start};
    $gffLine.=join("\t",$seqid,$source,$$s{type},$$s{start},$$s{end},$score,$strand,$phase,$$s{attributes})."\n";

    # update DR for the next iteration
    $DRstart=$$s{end}+1;

    # DC - reset $chaineDR
    $chaineDR = ""; # DC
  }
  # one more DR
  #$gffLine.=join("\t",$seqid,$source,"CRISPRdr",$DRstart,($DRstart+$tag{DR_length}-1),$score,$strand,$phase,"Parent=$featureId")."\n";
  
  $sequenceDR =  $dbSeq->seq($seqid, $DRstart => ($DRstart+$tag{DR_length}-1)); # DC - DR sequence
  $gffLine.=join("\t",$seqid,$source,"CRISPRdr",$DRstart,($DRstart+$tag{DR_length}-1),$score,$strand,$phase,"sequence=$sequenceDR;Parent=$featureId;ID=DR_$DRstart")."\n"; # DC - 05/2017 Adding DR sequence

  # DC - Right Flanking sequence
  $flankAttributes = "sequence=$rightFlankSeq;at%=".$rightAT; # DC
  $gffLine.=join("\t",$seqid,$source,"RightFLANK",($end+1),($end+$flankingRegion),$score,$strand,$phase,$flankAttributes.$leader2.";Parent=$featureId;ID=fl_$end")."\n"; # DC - show results

  $chaineDR = ">DR_".$DRstart."_".($DRstart+$tag{DR_length}-1); #DC
  print DRLIST "$chaineDR\n$sequenceDR\n"; #DC

  close DRLIST; # DC

  return ($gffLine,$idCrispr,$potentialDir,$globalAT);
}
#------------------------------------------------------------------------------
# DC - 05/2017 - makeJson (creates json file from GFF file)
sub makeJson
{
  my ($gff,$ResultDir,$refSeq) = @_; # retrieve gff file as returned by makeGff function
  my $input_spacer = $ResultDir."/".$refSeq."/Spacers_"; # Path to the Spacers files
  my $input_dr = $ResultDir."/".$refSeq."/DRs_"; # Path to the DRs files

  my $jsonCRISPRdata = ""; # variable representing all the data related to the former json file
  

  my ($crisprID,$name,$trusted,$score,$orientation,$potentiallyFalse,$start,$end,$DRconsensus,$DRlength,$nbSpacers,$sequence,$entropyDRs,$conservationSpacers,$leader,$leftAT,$rightAT,$meanSpacers,$ratioDRspacer);
  my $jsonLine = ""; # lines which will be written in JSON file
  my $secondChain = "";
  my ($file, $extension) = split(/\./, $gff); # retrieve file name
  # Write Json file
  #my ($hour,$min,$sec) = Now();

  if($logOption){
  	print LOG "\n Create JSON format data corresponding to CRISPR arrays\n"; #[$hour:$min:$sec]
  }
  
  # Hash consensusDR <=> evidence-level
  my %hashDRlevel = ();

  if($classifySmall){
  ## Pre-treatment in order to detect evidence-level 1 arrays having same consensus DR with an evidence-level 4 array.
  # Read and parse GFF file
  open(PRE,"<",$gff) or die("Could not open the GFF file $gff because $!\n");
  
  while(<PRE>) {
    chomp($_);

    if ($_ !~  m/^#/ && $_ =~  m/CRISPRCasFinder/) {
	
	my @tab = split(/\t/, $_);

	#### main CRISPR
        if($tab[2] eq "CRISPR" || $tab[2] eq "PossibleCRISPR"){

	  my %hash = ();

	  my @newTab = split(/;/, $tab[8]);
	  for (my $i = 0; $i <= $#newTab; $i++) {
	    my @otherTab = split(/=/, $newTab[$i]);
	    $hash{$otherTab[0]} = $otherTab[1];
	  }
	  #$DRconsensus = $hash{DR};
	  #$nbSpacers = $hash{Number_of_spacers};
	  $crisprID = $hash{ID};

	  ## Get Id
	  my @tabIdCRISPR = split(/_/, $crisprID);
	  my $idNumber = pop(@tabIdCRISPR);

	  ## Entropy DRs
	  if(-e $input_dr.$idNumber){
	  	my $drFasta = fastaAlignmentMuscle($input_dr.$idNumber);
		
	  	$entropyDRs = entropy($drFasta);
	  }

	  ## Conservation Spacers
	  if(-e $input_spacer.$idNumber) {
	  	my $spacerFasta = fastaAlignmentMuscle($input_spacer.$idNumber);
		$conservationSpacers = sequenceAlignmentMuscle($spacerFasta);
		$score .= "_".$conservationSpacers;

		my @tabSpacerLength=();
		open (SPACER, "<".$input_spacer.$idNumber) or die "open : $!";  #Open Spacer file

		while (defined (my $lineSp = <SPACER>)){
			chomp($lineSp);
			my $lengthInfo = 0;
			if($lineSp !~ /^>/){
				#$lineSp =~ s/\s//g;   #Remove spaces if exist
				$lengthInfo = length($lineSp);
				push(@tabSpacerLength,$lengthInfo);
			}
		}
		close SPACER;


	  }

	  if( ($entropyDRs >= 70) and ($conservationSpacers <= 8) and ($hash{Number_of_spacers} > 3) ){ 
		$hashDRlevel{$hash{DR}} = 4;
	  } 
        }
    }    
  }
  ## End of Pre-treatment
  close(PRE);
  }

  # Read and parse GFF file
  open(GFF,"<",$gff) or die("Could not open the GFF file $gff because $!\n");
  
  #open(JSON,">",$file."_crispr.json") or die("Could not open the JSON file $file _crispr.json because $!\n");
  #print JSON "[\n"; #replacement of "{" by "["
  $jsonCRISPRdata = "[\n";

  while(<GFF>) {
    chomp($_);
   
    if ($_ !~  m/^#/ && $_ =~  m/CRISPRCasFinder/) {
        
        my @tab = split(/\t/, $_);

	#### LeftFLANK (left flanking sequence)
	if($tab[2] eq "LeftFLANK"){
	  $start = $tab[3];
	  $end = $tab[4];
	  my @newTab = split(/;/, $tab[8]);
	  $sequence = $newTab[0];
	  $sequence =~ s/sequence=//;
	  $leftAT = $newTab[1];
   	  $leftAT =~ s/at%=//;
	  if($newTab[2] eq "leader"){
	    $leader = 1;
	  }
	  else {
	    $leader = 0;
	  }
	  $secondChain = "\"Regions\": [\n{\n\"Type\": \"LeftFLANK\",\"Start\": ".$start.",\"End\": ".$end.",\"Sequence\": \"".$sequence."\",\"Leader\": ".$leader.",\"AT\": ".$leftAT."\n},";
	}

	#### main CRISPR
        if($tab[2] eq "CRISPR" || $tab[2] eq "PossibleCRISPR"){

	  my %hash = ();
	  $start = $tab[3];
	  $end = $tab[4];

	  $orientation = $tab[6];
	  if($tab[6] eq "."){ $orientation = "ND"; }
	  if($tab[2] eq "CRISPR"){
	  	$trusted = 1;
	  }
	  else{
		$trusted = 0;
	  }
	  my @newTab = split(/;/, $tab[8]);
	  for (my $i = 0; $i <= $#newTab; $i++) {
	    my @otherTab = split(/=/, $newTab[$i]);
	    $hash{$otherTab[0]} = $otherTab[1];
	  }
	  $DRconsensus = $hash{DR};
	  $DRlength = $hash{DR_length};
	  $nbSpacers = $hash{Number_of_spacers};
	  $crisprID = $hash{ID};
	  $name = $hash{name};

	  ## Get Id
	  my @tabIdCRISPR = split(/_/, $crisprID);
	  my $idNumber = pop(@tabIdCRISPR);

	  ## Entropy DRs
	  if(-e $input_dr.$idNumber){
	  	my $drFasta = fastaAlignmentMuscle($input_dr.$idNumber);
		
	  	$entropyDRs = entropy($drFasta);
		#$score = "$entropyDRs";
	  }

	  ## Conservation Spacers
	  if(-e $input_spacer.$idNumber) {
	  	my $spacerFasta = fastaAlignmentMuscle($input_spacer.$idNumber);
		$conservationSpacers = sequenceAlignmentMuscle($spacerFasta);

		$score .= "_".$conservationSpacers;

		my @tabSpacerLength=();
		open (SPACER, "<".$input_spacer.$idNumber) or die "open : $!";  #Open Spacer file

		while (defined (my $lineSp = <SPACER>)){
			chomp($lineSp);
			my $lengthInfo = 0;
			if($lineSp !~ /^>/){
				#$lineSp =~ s/\s//g;   #Remove spaces if exist
				$lengthInfo = length($lineSp);
				push(@tabSpacerLength,$lengthInfo);
			}
		}
		close SPACER;

		##add Mean Spacers
		$meanSpacers = &average(\@tabSpacerLength);
		#$medianSpacers = median(@tabSpacerLength);
		##add ratio DR/Spacer
		$ratioDRspacer = $DRlength / $meanSpacers;


	  }
	  
	  ## Repeat ID + CRISPRDirection
	  my @tabR = repeatIDandNb($DRconsensus);
      	  my $crisprDirection = repeatDirection($tabR[0]);
	  if($crisprDirection =~  m/^F/){
		$crisprDirection = "+";
	  }
	  elsif($crisprDirection =~  m/^R/){
		$crisprDirection = "-";
	  }
	  else{
		$crisprDirection = "ND";
	  }

	  ## Evidence Level
	  my $eL = 0; # evidenceLevel (eL) = 4 means good CRISPR; eL = 3 means acceptable CRISPR; eL = 2 means bad CRISPR; eL = 1 means hypothetical CRISPR; 

	  if( ($entropyDRs >= 70) and ($conservationSpacers <= 8) and ($nbSpacers > 3) ){ $eL=4;}
	  if( ($entropyDRs >= 70) and ($conservationSpacers > 8) and ($nbSpacers > 3) ){ $eL=3;}
	  if( ($entropyDRs < 70) and ($nbSpacers > 3) ){ $eL=2;}
	  if( $nbSpacers <= 3 ){
		if(exists($hashDRlevel{$DRconsensus}) and $classifySmall){
			$eL=4;
			#print CRISPR into file "smallArraysReclassification.xls"
			my $smallArrays = $ResultDir."/smallArraysReclassification.xls";
			open (SMALLAR, ">>$smallArrays") or die "open : $!"; 
			print SMALLAR "$crisprID\t$DRconsensus\t1\t4\n";
			close (SMALLAR);
		}
		else{ 
	  		$eL=1;
		}
	  }

	  #if( ($entropyDRs >= 70) and ($conservationSpacers <= 8) and ($nbSpacers > 3) and ($ratioDRspacer>=0.8) and ($ratioDRspacer<=1.2) ){ $eL=4;}
      	  #if( ($entropyDRs >= 70) and ($nbSpacers > 3) and (($conservationSpacers > 8) or ($ratioDRspacer<0.8) or ($ratioDRspacer>1.2) ) ){ $eL=3;}
      	  #if( ($entropyDRs < 70) and ($nbSpacers > 3) ){ $eL=2;}
      	  #if( $nbSpacers <= 3 ){ $eL=1;}

	  ## Change $crisprID
	  $crisprID =~ s/PossibleCrispr_//; #PossibleCrispr_  Crispr_
	  $crisprID =~ s/Crispr_//;
	  
	  $jsonLine .= "{\n\"Name\": \"".$crisprID."\",\"Start\": ".$start.",\"End\": ".$end.",\"DR_Consensus\": \"".$DRconsensus."\",\"Repeat_ID\": \"".$tabR[0]."\",\"DR_Length\": ".$DRlength.",\"Spacers\": ".$nbSpacers.",\"Potential_Orientation\": \"".$orientation."\",\"CRISPRDirection\": \"".$crisprDirection."\",\"Evidence_Level\": ".$eL.",\"Conservation_DRs\": ".$entropyDRs.",\"Conservation_Spacers\": ".$conservationSpacers.",".$secondChain." " ;   # Please note that {\"value\": and } have been removed
	 #"}," 
	  
        }
	
	#### DR or Spacer
	if($tab[2] eq "CRISPRdr" || $tab[2] eq "CRISPRspacer"){
	  $start = $tab[3];
	  $end = $tab[4];
	  my @newTab = split(/;/, $tab[8]);
	  $sequence = $newTab[0];
	  $sequence =~ s/sequence=//;
	  my $seqType = "";
	  if($tab[2] eq "CRISPRdr"){
	    $seqType = "DR";
	  }
	  elsif($tab[2] eq "CRISPRspacer"){
	    $seqType = "Spacer";
	  }
	  
	  $secondChain = "{\n\"Type\": \"".$seqType."\",\"Start\": ".$start.",\"End\": ".$end.",\"Sequence\": \"".$sequence."\"\n},";
	  $jsonLine .= $secondChain;
	}

	#### RightFLANK
	if($tab[2] eq "RightFLANK"){
	  $start = $tab[3];
	  $end = $tab[4];
	  my @newTab = split(/;/, $tab[8]);
	  $sequence = $newTab[0];
	  $sequence =~ s/sequence=//;
	  $rightAT = $newTab[1];
   	  $rightAT =~ s/at%=//;
	  if($newTab[2] eq "leader"){
	    $leader = 1;
	  }
	  else {
	    $leader = 0;
	  }
	  $secondChain = "{\n\"Type\": \"RightFLANK\",\"Start\": ".$start.",\"End\": ".$end.",\"Sequence\": \"".$sequence."\",\"Leader\": ".$leader.",\"AT\": ".$rightAT."\n}\n]\n},";
	  $jsonLine .= $secondChain;
	}
	
    }
  }

  close GFF;

  $jsonLine .= "]\n"; # replacement of "}" by "]"
  $jsonLine =~ s/,]/]/; # change end of file to better fit JSON format
  
  $jsonLine =~ s/,/,\n/g; # make file more readable
  #print JSON $jsonLine;

  $jsonCRISPRdata .= $jsonLine;

  #close JSON;
  return $jsonCRISPRdata;
  
}
#------------------------------------------------------------------------------
# DC - 06/2017 - foundInCRISPRdb
sub foundInCRISPRdb
{
  my ($seq,$start,$end) = @_;
  my $found = 0; # boolean value to indicate if the same CRISPR has been found in DB.
  my %hashSeq = ();
  #my %hashStart = ();
  #my %hashEnd = ();
  my ($dbSeq,$dbStart,$dbEnd,$line); 

  #Translate $seq
  my @seqID1 = (); # DC
  if ($seq =~ /\|/){
    @seqID1 = split (/\|/,$seq);
    $seq = pop(@seqID1);
    @seqID1 = ();
    @seqID1 = split (/\./,$seq);
    $seq = $seqID1[0];
  }
  else{
    if($seq =~ /\./){
      @seqID1 = split (/\./,$seq);
      $seq = $seqID1[0];
    }
  }  

  #my $currentRepository = getcwd();
  my $file = $crisprdb;  #"all_CRISPRs_DB.csv"; # file containing alidated CRISPRs from DB

  if (-e $file){ 
	  open (FILE, "<$file") or die "open $file: $!";  #Open $file
	  #print "SEQ=$seq ; Start=$start ; End=$end\n";

	  while (defined ($line = <FILE>)){
		chomp($line);
		#$line =~ s/\s//g;
	    if($line !~  m/^Specie\tRefSeq/){ # DC - 07/2017 - check if file line is not beginning by "Specie\tRefSeq"
		my @tab = split(/\t/,$line);
	
		($dbSeq,$dbStart,$dbEnd) = ($tab[1],$tab[4],$tab[5]);
		#print "DBStart = $dbStart ; DBEnd = $dbEnd\n";
		my $key = $dbSeq."_".$dbStart."_".$dbEnd;
		my $secondKey = $dbSeq."_".$dbStart."_".($dbEnd+1); #in case the end of sequence is $dbEnd+1
		my $thirdKey = $dbSeq."_".$dbStart."_".($dbEnd-1); #in case the end of sequence is $dbEnd-1

		$hashSeq{$key}++;
		$hashSeq{$secondKey}++;
		$hashSeq{$thirdKey}++;
	    }
	  }
	  close (FILE);
  }

  my $key2 = $seq."_".$start."_".$end;
  #print "Start = $start ; End = $end\n";
  if($hashSeq{$key2}){
    $found = 1
  }
  else{
    $found = 0;
  }

  return $found;
}
#------------------------------------------------------------------------------
# DC - 05/2017 - repeatID (CRISPRdb) + Expert based orientation
sub repeatIDandNb
{
  my ($repeat) = @_;
  chomp($repeat);
  $repeat =~ s/\s//g; # replace all spaces by nothing

  my $file = $repeats; #"Repeat_List.csv"; # file containing repeat sequences, IDs, and Nb in CRISPRdb (last update)
  my $line = "";
  my ($id,$number,$temp,$orientation); #repeatID and repeat number occurred in CRISPRdb, +orientation
  my %hashID = ();
  my %hashNB = ();
  my %hashOrientation = ();

  if (-e $file){
	  open (FILE, "<$file") or die "open $file: $!";  #Open $file

	  while (defined ($line = <FILE>)){
		chomp($line);
		$line =~ s/\s//g;
		($temp,$id,$number,$orientation) = split(/;/,$line);

		$hashID{$temp} = $id;
		$hashNB{$temp} = $number;
		$hashOrientation{$temp} = $orientation;
	  }
	  close (FILE);
  }

  if ($hashID{$repeat}) {
	$id = $hashID{$repeat};
	$number = $hashNB{$repeat};
	$orientation = $hashOrientation{$repeat};
  }
  else {
	$id = "Unknown";
	$number = 0;
	$orientation = "";
  }
  return ($id,$number,$orientation);
}
#------------------------------------------------------------------------------
# DC - 05/2017 - retrieve results from CRISPRDirection tool 
sub repeatDirection
{
  my ($id) = @_;
  chomp($id);
  $id =~ s/\s//g; # replace all spaces by nothing

  my $file = $dirRepeat; #"repeatDirection.txt"; # file containing repeat IDs and Orientation according to CRISPRDirection
  my $line = "";
  my ($orientation,$temp); #orientation and id
  my %hashID = ();

  if (-e $file){
	  open (FILE, "<$file") or die "open $file: $!";  #Open $file
	  while (defined ($line = <FILE>)){
		chomp($line);
		($temp,$orientation) = split(/\t/,$line);
		chomp($temp);
		$temp =~ s/\s//g;
		$hashID{$temp} = $orientation;
	  }
	  close (FILE);
  }
  if ($hashID{$id}) {
	return $hashID{$id};
  }
  else {
	return "Unknown";
  }
}
#------------------------------------------------------------------------------
# DC - 05/2017 - analyzeCRISPR (provide some stats)
sub crisprAnalysis
{
  my($ResultDir,$refSeq,$nbCrispr,$seqDescription,@idDir)=@_;

  my ($seqid,$start,$end,$potentiallyFalse,$DRconsensus,$DRlength,$nbSpacers);
  my $input_file_confirmed = $ResultDir."/".$refSeq."/".$refSeq."_Crispr_"; # Path to confirmed CRISPRs files
  my $input_spacer = $ResultDir."/".$refSeq."/Spacers_"; # Path to the Spacers files
  my $input_dr = $ResultDir."/".$refSeq."/DRs_"; # Path to the DRs files
  my $input_file_possible = $ResultDir."/".$refSeq."/".$refSeq."_PossibleCrispr_"; # Path to possible/hypothetical CRISPRs files
  my %hashIdDir = ();

  #Add chdir current repository DC
  #chdir($currentRepository); #DC
	
  my $directory = $ResultDir."/".$refSeq."/Alignments"; # directory  to store Alignments
  unless(-d $directory){ mkdir $directory or die "$0: I can not create the folder $directory: $!\n" }

  ## Test if Seq description/name is defined
  if (defined $seqDescription && $seqDescription ne '') {
    # do nothing
  }
  else {
    $seqDescription = "NA";
  }

  ### Hash Id_Dir (ID of crisprs and their potential direction)
  if(@idDir) {
    foreach my $v (@idDir) {
	chomp ($v);
	my ($id,$direction) = split(/_/,$v);  
    	$hashIdDir{$id} = $direction;  
    }
  }
  ### Write Results
  #my ($hour,$min,$sec) = Now();

  if($logOption){
  	print LOG "\n Write Crisprs_REPORT.tsv\n"; #[$hour:$min:$sec]
  }

  my $resultsCRISPRs = $ResultDir."/Crisprs_REPORT.tsv";
  my $resultsCRISPRCasSummary = $ResultDir."/CRISPR-Cas_summary.tsv";  # TSV file which will contain a summary for CRISPRs and Cas
  my @tabCCScrispr = (); # table to store Crisprs


  my $header = "Strain\tSequence\tSequence_basename\tDuplicated_Spacers\tCRISPR_Id\tCRISPR_Start\tCRISPR_End\tCRISPR_Length\tPotential_Orientation (AT%)\tCRISPRDirection\tConsensus_Repeat\tRepeat_ID (CRISPRdb)\tNb_CRISPRs_with_same_Repeat (CRISPRdb)\tRepeat_Length\tSpacers_Nb\tMean_size_Spacers\tStandard_Deviation_Spacers\tNb_Repeats_matching_Consensus\tRatio_Repeats_match/TotalRepeat\tConservation_Repeats (% identity)\tEBcons_Repeats\tConservation_Spacers (% identity)\tEBcons_Spacers\tRepeat_Length_plus_mean_size_Spacers\tRatio_Repeat/mean_Spacers_Length\tCRISPR_found_in_DB (if sequence IDs are similar)\tEvidence_Level\n";
  #adding stuff after Conservation_Spacers (% identity): i.e EBcons_Spacers, Repeat_Length_plus_mean_size_Spacers, ratio_Repeat_Length_mean_size_Spacers.

  my $totalFound = 0;
  my $alreadyExist = 0;

  my %hashEvidence =(); # Hash table for counting evidence-levels
  $hashEvidence{1}=0;
  $hashEvidence{2}=0;
  $hashEvidence{3}=0;
  $hashEvidence{4}=0;

  if(-e $resultsCRISPRs){
    $alreadyExist = 1;
  }
  else{
    $alreadyExist = 0;
  }

  open (RESULTCCS, ">>$resultsCRISPRCasSummary") or die "open : $!";  #Open $resultsCRISPRCasSummary
  print RESULTCCS "$refSeq\t";

  open (RESULT, ">>$resultsCRISPRs") or die "open : $!";  #Open $resultsCRISPRs

  ## Write the title line in this file (if it not already exists)
  
  if ($alreadyExist == 0) {
    print RESULT $header;
  }
  else {
    print RESULT "\n";
  }

  ## Analyze by CRISPR file
  # $j is a 'counting' variable for (possible/hypothetical or confirmed) CRISPR file

  for (my $j = 1; $j <= $nbCrispr; $j++) {

    my $spacerFasta = "";
    my $drFasta = "";
    my $crispr = "";
    my $hypothetical = 0;
    
    if(-e $input_file_confirmed.$j){
      ($seqid,$start,$end,$potentiallyFalse,$DRconsensus,$DRlength,$nbSpacers) = analyzeCrisprs ($input_file_confirmed.$j);
      $crispr = "CRISPR";

      if($nbSpacers <= 3){
        $hypothetical = 1;
      }
    }
    elsif (-e $input_file_possible.$j) {
      ($seqid,$start,$end,$potentiallyFalse,$DRconsensus,$DRlength,$nbSpacers) = analyzeCrisprs ($input_file_possible.$j);
      $crispr = "PossibleCRISPR";

      if($nbSpacers <= 3){
        $hypothetical = 1;
      }
    }
    
    ##Found in DB?
    my $found = foundInCRISPRdb($seqid,$start,$end);
    if($found){
   	$totalFound++;     
    }
        #Manage Mean and Standard Deviation of spacers size And Conversion
	##read fasta file containing spacers
	## ($spacerLength, $meanSpacers, $standardSpacers)
	my (@tabSpacerLength, $meanSpacers, $medianSpacers, $standardSpacers, $identity, $entropy, $identityDR, $entropySpacers, $DRspacer, $ratioDRspacer);     

     #### Analyze DRs
        ## Verify how many DRs actually match with the consensus DR sequence
	my $countDRmatch =0;    # to count DRs matching with consensus
	my $countAllLinesDR =0;   #To count all lines in DR fasta file

	if (-e $input_dr.$j){

		my $lineDR="";
		open (DRFILE, "<$input_dr".$j) or die "open : $!";  #Open DR file

		while (defined ($lineDR = <DRFILE>)){
			chomp($lineDR);
			$countAllLinesDR++;
			if($lineDR !~ /^>/){
				if($lineDR =~ /$DRconsensus/){
					$countDRmatch++;
				}
			}
		}
		close (DRFILE);

		## Entropy DRs
	
			$drFasta = fastaAlignmentMuscle($input_dr.$j);
			$entropy = entropy($drFasta);
			$identityDR = sequenceAlignmentMuscle($drFasta);
	}
    #### End Analyze DRs

    #### Analyze Spacers

	if (-e $input_spacer.$j){

		my $lineSp="";
		my %duplicatedSpacer=();
		open (SPACER, "<".$input_spacer.$j) or die "open : $!";  #Open Spacer file

		while (defined ($lineSp = <SPACER>)){
			chomp($lineSp);
			my $lengthInfo = 0;
			if($lineSp !~ /^>/){
				#$lineSp =~ s/\s//g;   #Remove spaces if exist
				$lengthInfo = length($lineSp);
				push(@tabSpacerLength,$lengthInfo);
				$duplicatedSpacer{$lineSp}++;
			}
		}
		close SPACER;
		
		##Check for repeated spacers
		foreach my $k2 (keys(%duplicatedSpacer)) {
			if($duplicatedSpacer{$k2} > 1){
				$potentiallyFalse = 1;
			}
		}
		##

		##add Mean and Standard deviation
		#print "@tabSpacerLength\n";
		$meanSpacers = &average(\@tabSpacerLength);
		$standardSpacers = &stdev(\@tabSpacerLength);
	 	#$medianSpacers = median(@tabSpacerLength);
		##add DR+Spacer length and ratio
		$DRspacer = $DRlength + $meanSpacers;
		$ratioDRspacer = $DRlength / $meanSpacers;

		#### %Alignment spacers
		
			$spacerFasta = fastaAlignmentMuscle($input_spacer.$j);
			$identity = sequenceAlignmentMuscle($spacerFasta);

		#Entropy based conservation of spacers

		#$drFasta = fastaAlignmentMuscle($input_dr.$j);
		$entropySpacers = entropy($spacerFasta);
	}
    #### End Analyze Spacers

    #### Move created alignment files into $directory
      #my ($hourMv,$minMv,$secMv) = Now();
      my $fastaStar = $spacerFasta."*";
      my $ret = system("cp -r $fastaStar $directory");
      #print "Moving files $fastaStar in $directory (return type: $ret)\n";
      
      if($logOption){	
      	print LOG "\n Moving files $fastaStar in $directory (return type: $ret)\n"; #[$hourMv:$minMv:$secMv]
      }

      $fastaStar = $drFasta."*";
      $ret = system("cp -r $fastaStar $directory");
      #print "Moving files $fastaStar in $directory (return type: $ret)\n";

      if($logOption){
      	print LOG "\n Moving files $fastaStar in $directory (return type: $ret)\n"; #[$hourMv:$minMv:$secMv]
      }

    #### Write results    
      my @tabR = repeatIDandNb($DRconsensus);
      my $crisprDirection = repeatDirection($tabR[0]);
      my $crisprID = $refSeq."_".$j; 	# DC removed ."_".$crispr 
      my $crisprLength = $end - $start;
      my $ratioDrsMatch = $countDRmatch / ($nbSpacers+1);
      #my $potentiallyTrue = 0; # probable decision (if good entropyDRs or good conservationSpacers, we could admit that the CRISPR is good, i.e it is = 1)

      my $eL = 0; # evidenceLevel (eL) = 4 means good CRISPR; eL = 3 means acceptable CRISPR; eL = 2 means bad CRISPR; eL = 1 means hypothetical CRISPR; 

      if( ($entropy >= 70) and ($identity <= 8) and ($nbSpacers > 3) ){ $eL=4;}
      if( ($entropy >= 70) and ($nbSpacers > 3) and ($identity > 8) ){ $eL=3;}
      if( ($entropy < 70) and ($nbSpacers > 3) ){ $eL=2;}
      if( $nbSpacers <= 3 ){ $eL=1;}

      #if( ($entropy >= 70) and ($identity <= 8) and ($nbSpacers > 3) and ($ratioDRspacer>=0.8) and ($ratioDRspacer<=1.2) ){ $eL=4;}
      #if( ($entropy >= 70) and ($nbSpacers > 3) and (($identity > 8) or ($ratioDRspacer<0.8) or ($ratioDRspacer>1.2) ) ){ $eL=3;}
      #if( ($entropy < 70) and ($nbSpacers > 3) ){ $eL=2;}
      #if( $nbSpacers <= 3 ){ $eL=1;}


      if ($eL >= $levelMin){
      print RESULT "$seqDescription\t$seqid\t$refSeq\t$potentiallyFalse\t$crisprID\t$start\t$end\t$crisprLength\t$hashIdDir{$j}\t$crisprDirection\t";
      print RESULT "$DRconsensus\t$tabR[0]\t$tabR[1]\t$DRlength\t$nbSpacers\t$meanSpacers\t$standardSpacers\t$countDRmatch\t$ratioDrsMatch\t";
      print RESULT "$identityDR\t$entropy\t$identity\t$entropySpacers\t$DRspacer\t$ratioDRspacer\t$found\t$eL\n";
      }
      #print in RESULTCCS and push CRISPR info in table @tabCCScrispr
      
      my $ccsCRISPR = $crisprID."[".$start.";".$end."] "."(evidence-level=".$eL."),";
      push (@tabCCScrispr, $ccsCRISPR);
      
      #Hash table containing counts of each evidence-level for the given sequences
      $hashEvidence{$eL}++;
      
      ## Write raw CRISPRs sequences
       # get subsequence
       my $inputfileSeq = $refSeq.".fna";
       if(-e $inputfileSeq){
		my $crisprSeqFile = $ResultDir."/rawCRISPRs.fna";
		open (RAWCRISPRS, ">>$crisprSeqFile") or die "open : $!";  #Open crisprSeqFile

	      	my $dbSeq = Bio::DB::Fasta->new($inputfileSeq); # sequence file (or $userfile)
	       
	       	# DC - 11/2017 - get raw sequence of CRISPR
	       	my $rawCRISPRsequence = $dbSeq->seq($seqid, $start => $end); # CRISPR sequence

		my $len = 80;

		my $formatted_seq = ">".$crisprID." ".$start.",".$end."\n";
		while (my $chunk = substr($rawCRISPRsequence, 0, $len, "")) {
		        $formatted_seq .= "$chunk\n";
		}

		print RAWCRISPRS $formatted_seq; # ."\n";

		close(RAWCRISPRS);
       }
      ##

  }

  #if($crisprdb eq ""){
  	#print RESULT "(Total_Nb_CRISPRs = $nbCrispr)\n";
  #}
  #elsif(-e $crisprdb){
  #	print RESULT "(Nb_CRISPRs_found_in_CRISPRdb/Total_CRISPRs = $totalFound / $nbCrispr)\n";
  #}
  
  print RESULT "\n";

  close RESULT;

  # string of evidence-levels count
  #foreach my $keyEvidence (keys %hashEvidence)
  #{
  	# if value is not defined put 0
  #	if ($hashEvidence{$keyEvidence} eq "") { $hashEvidence{$keyEvidence} = 0; }
  #}

  my $evidenceString = "Nb_arrays_evidence-level_1=".$hashEvidence{1}.",Nb_arrays_evidence-level_2=".$hashEvidence{2};
  $evidenceString .= ",Nb_arrays_evidence-level_3=".$hashEvidence{3}.",Nb_arrays_evidence-level_4=".$hashEvidence{4};

  # print table of CRISPRs in RESULTCCS
  if($launchCasFinder){
    if($rcfowce and ($nbCrispr==0)){
      print RESULTCCS "@tabCCScrispr\t$nbCrispr\tND\n"; # addition of evidence-level info
    }
    else{
      print RESULTCCS "@tabCCScrispr\t$nbCrispr\t$evidenceString\t"; # addition of evidence-level info
    }
  }
  else{
    print RESULTCCS "@tabCCScrispr\t$nbCrispr\t$evidenceString\n"; # addition of evidence-level info
  }


  close RESULTCCS;

  return ($resultsCRISPRs,$totalFound);
}
#------------------------------------------------------------------------------
sub analyzeCrisprs
{
  my($report)=@_;
  my($seqid,$start,$end,$strand,$attributes,%tag,$DRconsensus,$DRlength,$nbSpacers);
  my $potentiallyFalse = 0;

  open(REPORT,"<",$report) or die("Could not open the report file $report because $!\n");

  while(<REPORT>){
    my($spacerStart,$spacerEnd,$spacerSequence);

    if(/#Potentially false/){ # DC - control if the given CRISPR is potentially false (#Potentially )
      $potentiallyFalse = 1; # boolean value to indicate false CRISPR
    }

    if(/Id: (\S+)/){
      $seqid=$1;
    }
    elsif (/Crispr.+?(\d+)\s+Crispr\D+(\d+)/){
      $start=$1;
      $end=$2;
      #print "START = $start; END = $end\n";
      my $nextLine=<REPORT>;
      $nextLine=~s/://g; # remove colons
      my @field=split(/\s+/,$nextLine);
      shift(@field); # remove the hash sign
      for(my $i=0;$i<@field;$i+=2){
        $attributes.="$field[$i]=$field[$i+1];";
        $tag{$field[$i]}=$field[$i+1];
  	#DC
	#print "ATTRIBUTES = $attributes\n";
      }
    }
    
  }
  close REPORT;
  
  if($tag{'DR'}) { $DRconsensus = $tag{'DR'}; }
  else { $DRconsensus = "NA"; }

  if($tag{'DR_length'}) { $DRlength = $tag{'DR_length'}; }
  else { $DRlength = "NA"; }

  if($tag{'Number_of_spacers'}) { $nbSpacers = $tag{'Number_of_spacers'}; }
  else { $nbSpacers = "NA"; }

  return ($seqid,$start,$end,$potentiallyFalse,$DRconsensus,$DRlength,$nbSpacers);
}
#------------------------------------------------------------------------------
sub compare
{
  my ($val1, $val2) = @_;
  return ($val1<=$val2*0.5) || ($val1>=$val2 * 1.5) ? 1 : 0;
  #if( ($val1<=$val2*0.5) || ($val1>=$val2 * 1.5) ){return 1;} else {return 0;}
}
#------------------------------------------------------------------------------
# DC - 05/2017 - Calculate AT%
sub atpercent
{
  my ($seq) = @_;
  my @charSeq = split(//, uc($seq));
  
  my %hashFlank = ();

  foreach my $v (@charSeq) {
    $hashFlank{$v} += 1;  
  }
  
  if(! $hashFlank{'A'}){ $hashFlank{'A'} = 0; }
  if(! $hashFlank{'T'}){ $hashFlank{'T'} = 0; }

  if(length($seq) == 0){
	return 0;
  }
  else{
  	return (($hashFlank{'A'} + $hashFlank{'T'})/(length($seq))) * 100;
  }
}

#------------------------------------------------------------------------------
sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

#------------------------------------------------------------------------------
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

#------------------------------------------------------------------------------
sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub max {
    splice(@_, ($_[0] > $_[1]) ? 1 : 0, 1);
    return ($#_ == 0) ? $_[0] : max(@_);
}

sub min {
    splice(@_, ($_[0] > $_[1]) ? 0 : 1, 1);
    return ($#_ == 0) ? $_[0] : min(@_);
}
#------------------------------------------------------------------------------
#DC - get %conservation (overall % identity) from a Muscle sequence alignment 
sub sequenceAlignmentMuscle
{
  my($file) = @_;
  my $ident = -1;
  eval
  {
	my $str = Bio::AlignIO->new(-file => $file);
  	my $aln = $str->next_aln();

	$ident = $aln->overall_percentage_identity; #overall_percentage_identity;
	
  };
 
  if($@)
  {
     	# An error occurred...
	$ident = -1;
  }

  return $ident;
}

# DC - 05/2017 - entropy calculation
#------------------------------------------------------------------------------
sub log2 {
  my $n = shift;
  return (log($n)/ log(2));
}

#------------------------------------------------------------------------------
# Function allowing to calculate conservation of DRs based on entropy
sub entropy
{
  my($file) = @_; #file given as parameter is an alignment Fasta file
  
  open F, $file or die "an error occurred while opening $file\n";
  my @lines = <F>;
  close F;

  my (%words, $total, @text);

  my $seqLength=0;
  my @tableRows=();
  my @tableCols=();

  my @table=();

  my $countLine = 0;

  foreach my $line (@lines) {
	chomp $line;
	my @words = split /[^a-zA-Z]+/, $line;
	#Replace to treat Fasta sequence
	#Retrieve each letter by column
	if ($line !~ />/)
	{
		$seqLength=length($line);
		#print "Line is : $line\n";
		for (my $i = 0; $i < $seqLength; $i++) {
			my $value = substr($line, $i, 1);
			$table[$countLine][$i] = $value;

			#Instanciate $tableCols with values from columns
			$tableCols[$i][$countLine] = $value;
		}		
		$countLine++;
	}
  }

  # entropy
  my $sumEntropy =0;

  for (my $j = 0; $j < $seqLength; $j++) {	
	my $sizeTable = @{$tableCols[$j]};
	my %element =();

	foreach my $elem (@{$tableCols[$j]}) {
		if ($elem eq "-") { 
		  $element{$elem} = 0; 
		} 
             	else{
		  $element{$elem}++;
		}
	}
	my $entropy;
	foreach my $word (@{$tableCols[$j]}) {
		if ($word ne "-"){
		  my $prob = $element{$word} / $sizeTable;
		  $entropy += log2($prob); 
		}
	}

	$entropy *= -1;
	$entropy = $entropy/$sizeTable;

	#Calculate Sum of entropies and divide by $seqLength ... Adding (1- $entropy) to get percentage-like results
	$sumEntropy = $sumEntropy + (1 - $entropy); 
  }

  #PRINT FINAL Entropy (Sum + "Mean")
  my $finalResult = ($sumEntropy / $seqLength) * 100; 
  
  if($finalResult<0){ $finalResult = 0; }
  
  return $finalResult;
}
#------------------------------------------------------------------------------
# Test if Prog is installed
sub isProgInstalled {
 
	my $program = shift;
 	my @PATH=split(":","$ENV{'PATH'}");     #Get users Path
	my $found = 0;
	foreach my $path (@PATH) 
	{
   	  if (-x "$path/$program") 
   	  {
   		  $found= 1;
   		  last;
    	  }
	}
	
	if ( $program =~ /^muscle/ )
	{
        	unless ($found) {
                	print "\nThe program $program cannot be found on your system\n";
			print "Have you installed it? Your PATH variable contains: $ENV{'PATH'}\n\n";
			print "\nPlease install $program\n";
			exit(1);
		}
 	}	
     return $found;
}

#------------------------------------------------------------------------------
sub extractcrispr
{
  my($seqfile,@spacers) = @_;
  my($count,$crisprfile,$seq,$spacerseq, $beg, $end, $id);

	$crisprfile = "spacers".$seqfile;
	open WRITER,"> $crisprfile" or die "The file cannot be edited !\n";
	my $in = new Bio::SeqIO(-format => 'Fasta',-file => $seqfile);
	$seq = $in->next_seq;

	$count=0;
	while($count<= $#spacers)
	{
		$beg = $spacers[$count]->Pos1;
		$end = $beg + $spacers[$count]->Length-1;
		$spacerseq = $seq->subseq($beg, $end);
		my $u = $count+1;
		$id = "spacer".$u;
		print WRITER ">".$id."\n";
		print WRITER $spacerseq."\n";
		$count++;
	}
	close WRITER;
}
#------------------------------------------------------------------------------
sub definespacers
{  #function that returns a structure of hypothetical crisprs definition (the structure spacers )
  my($DRlength, $crisprfile,$indexname)= @_;
  my(@lines, $count, @temp,$posdeb, $posend, $len,$i,$refFalsSpacers, @spacers_correct);
  my @spacers = ();
  my $nbspacers = 0;

my $mism = 0;
my $Toterr = 0;
my $simDRs = 0; # indicates if all DRs are different
  open(FD, $crisprfile)  or die "Error in opening the input file : $crisprfile!";
  $i = 0;	#counter for the number of spacers we find
  my $repCount = 0;
  my $repCurrent = 0;
  while(my $line = <FD>)
  {
     if($line=~/^# HitCount:\s+(\d+)/i){
        $repCount = $1+0;
     }
     next if $line =~ /^#/;
     $line = trim($line);
     next unless $line;
     next if $line =~ /^Start/;
     $repCurrent++;
     @temp = split(/ +/, $line);
     # store the first DR start end position
     unless(defined $posend){
        $posend = $temp[1];
	next;
     }
     $posdeb = $temp[0];
     $len = $posdeb - $posend -1 ;
     if($len < ($DRlength*$Sp1) ) # DC - replaced "if($len < ($DRlength*0.6) )" by "if($len < ($DRlength*$Sp1) )"
     {
	$mism = $temp[$#temp-1] eq "." ? 0 : $temp[$#temp-1]/$DRlength;
     }
     else
     {
	if( $repCurrent == $repCount )
	{
		if($len > ($DRlength*$Sp2)) # DC- replaced "if($len > ($DRlength*2.5))" by "if($len > ($DRlength*$Sp2))"
		{
			$posend = $temp[1];
		} 
		else
		{
			@spacers = add_spacer($posend,$len,$i,@spacers);
			$simDRs = 1 if $temp[$#temp-1] eq '.';
			$Toterr = $Toterr + $mism;
			$posend = $temp[1];
			$i++;
		}
		$mism = $temp[$#temp-1] eq "." ? 0 : $temp[$#temp-1]/$DRlength;
	}
	else
	{
		@spacers = add_spacer($posend,$len,$i,@spacers);
		$simDRs = 1 if $temp[$#temp-1] eq '.';
		$Toterr = $Toterr + $mism;
		$posend = $temp[1];
		$mism = $temp[$#temp-1] eq "." ? 0 : $temp[$#temp-1]/$DRlength;
		$i++;
	}
     }
   }
   close FD;
   
   $Toterr = ($Toterr + $mism)/($i+1);
   if($Toterr > $DRerrors) # DC - replaced if($Toterr > 0.2) by if($Toterr > $DRerrors)
   {
   	@spacers = ();
	my @false_spacers = ();
	$refFalsSpacers = \@false_spacers;
   }
   
   if($#spacers >=0)
   {
	extractcrispr($indexname,@spacers);
	($refFalsSpacers, @spacers_correct) = check_cris_div($DRlength,$indexname, @spacers);
   }
   else 
   {
   	my @false_spacers = ();
   	$refFalsSpacers = \@false_spacers;
   }

   #my ($hour,$min,$sec) = Now();

   if($logOption){
	print LOG " Function definespacers returns following values for ($crisprfile): $simDRs , $refFalsSpacers , @spacers_correct\n "; #[$hour:$min:$sec]
   }
   
   return($simDRs,$refFalsSpacers, @spacers_correct);
}

#------------------------------------------------------------------------------
# this function checks the spacers length (comparison with $Sp2 DR)
# deletes the wrong spacers in the begining 
# returns an array of wrong spacers begin positions 
sub check_cris_div
{
  my($DRlength,$indexname, @spacers) = @_;
  my @false_spacers = ();
  my($i,$j,$comp,$begtest,$spLen,$spPos, @spacers_new);
  @spacers_new = ();
  $j = 0;
  $begtest = 1;
  $comp = $DRlength * $Sp2; # DC - replaced "2.5" by "$Sp2"
  my $reffal_spacers;

  for($i=0; $i<= $#spacers; $i++)
  {
		$spLen = $spacers[$i]-> Length;
		if( ($spLen > $comp) )
		{# DC- uncommenting the whole IF condition

	  		if($i == $#spacers)
	  		{
				if($#spacers > 0)
				{
					my $crisprfile = "spacers".$indexname;
					open(FD, $crisprfile)  or die "Error in opening the input file : $crisprfile!";
					my @lines=<FD>;
					close(FD);
					open WRITER,"> TmpSpacers" or die "The file TmpSpacers cannot be edited !\n";
					my $k = 0;
					while($k<$#lines -1)
					{
		  				print WRITER $lines[$k];
		  				$k++;
					}
					close WRITER;
					rename("TmpSpacers",$crisprfile);
					# not good : @spacers_new is not the solution!!
					($reffal_spacers,@spacers_new)=check_cris_div($DRlength,$indexname, @spacers_new);
				}
	  		}
	  		else
	  		{
	    			if(!$begtest)
	    			{
					push(@false_spacers,$spacers[$i]-> Pos1);
					$spPos = $spacers[$i]-> Pos1 - 1;
					@spacers_new = add_spacer($spPos,$spLen,$j,@spacers_new);
					$j++; $begtest = 0;
				}
	    			else
	    			{
					my $crisprfile = "spacers".$indexname;
					open(FD, $crisprfile)  or die "Error in opening the input file : $crisprfile!";
					my @lines=<FD>;
					close(FD);
					open WRITER,"> TmpSpacers" or die "The file TmpSpacers cannot be edited !\n";
					my $k = 2;
					while($k<=$#lines)
					{
		  				print WRITER $lines[$k];
		  				$k++;
					}
					close WRITER;
					rename("TmpSpacers",$crisprfile);
	    			}
			}
		}
		elsif($spLen <= $comp)
		{
			$spPos = $spacers[$i]-> Pos1 - 1;
			@spacers_new = add_spacer($spPos,$spLen,$j,@spacers_new);
			$j++; $begtest = 0;
	   	}
   }
	if(!$reffal_spacers){$reffal_spacers = \@false_spacers;}

	#my ($hour,$min,$sec) = Now();

   	if($logOption){
		print LOG " Function check_cris_div returns: $reffal_spacers , @spacers_new\n "; #[$hour:$min:$sec]
   	}

	return($reffal_spacers,@spacers_new);
}
#------------------------------------------------------------------------------
sub add_spacer
{
  my($pos,$len,$i,@spacers)=@_;
  $pos = $pos+1;
  $spacers[$i] = Rep->new();
  $spacers[$i]-> Pos1($pos);
  $spacers[$i]-> Length($len);

  return @spacers;
}
#------------------------------------------------------------------------------

# sub checkspacersAlignMuscle (do the same as above but using Muscle)
sub checkspacersAlignMuscle
{
    my($file) = @_;
    my $ident;
    eval
    {
	my $resultAlign = fastaAlignmentMuscleOne($file);
      	
	my $str = Bio::AlignIO->new(-file => $resultAlign);
  	my $aln = $str->next_aln();

	$ident = $aln->percentage_identity;
    };
 
    if($@)
    {
        # An error occurred...
	$ident = 100;
    }

    if($ident >= $SpSim) {return 0;} else {return 1;} # DC - replaced 60 by $SpSim

}

#------------------------------------------------------------------------------
# Muscle alignment
sub fastaAlignmentMuscle
{
  my($file) = @_;
  my $result = $file."_fasta";

  eval
  {
	my $prog = isProgInstalled("muscle");
	my $muscle = " muscle -align $file -output $result ";
	#my $muscle = " muscle -in $file -out $result "; #muscle command-line
	if ($quiet){
		$muscle .= " -quiet ";
	}
	#if ($fast){
	#	$muscle .= " -maxiters 1 -diags ";
	#}
	
	#my ($hour,$min,$sec) = Now();

	if($prog){
		makesystemcall($muscle);

		if($logOption){
			print LOG "\n $muscle\n"; #[$hour:$min:$sec]
		}
	}
  };

  if($@)
  {
     	# An error occurred...
	print "An error occurred in function fastaAlignmentMuscle\n";
  }

  return $result;
}

#------------------------------------------------------------------------------
# Muscle alignment for specified purpose
sub fastaAlignmentMuscleOne
{
  my($file) = @_;
  my $result = "fastaMuscle_".$file;

  eval
  {
	my $prog = isProgInstalled("muscle");
	my $muscle = " muscle -align $file -output $result ";
	#my $muscle = " muscle -in $file -out $result "; #muscle command-line
	if ($quiet){
		$muscle .= " -quiet ";
	}
	#if ($fast){
	#	$muscle .= " -maxiters 1 -diags ";
	#}

	#my ($hour,$min,$sec) = Now();

	if($prog){
		makesystemcall($muscle);
	}
  };

  if($@)
  {
     	# An error occurred...
	print "An error occurred in function fastaAlignmentMuscle\n";
  }

  return $result;
}

#------------------------------------------------------------------------------

sub extractsequence
{
  # extract a genomic sequence from a well defined position till a predefined length
  #extractsequence(FILENAME, startPOS1, endPOS1, startPOS2, endPOS2,...)
  # the output is files containig each hypothetical crispr loci
  my(@seqencesdefinition) = @_;
  my($count,$seqfile, $number,@subselectoptions); 
  for($count = 0; $count < $#seqencesdefinition; $count+=2)
  {
		# extending the sequences delimiters and testing that they doesn't overflow
		#positions cases
    	if($seqencesdefinition[$count] - 500 > 0){$seqencesdefinition[$count]=$seqencesdefinition[$count] - 500;} 
    	else{$seqencesdefinition[$count]=0;}
		#if($seqencesdefinition[$count+1]+500> $seq->length()){$seqencesdefinition[$count+1]=$seq->length()-1;} # NV
		if($seqencesdefinition[$count+1]+500>= $seq->length()){$seqencesdefinition[$count+1]=$seq->length()-1;}
		else{ $seqencesdefinition[$count+1] = $seqencesdefinition[$count+1] + 500;}
    	$number = ($count+2)/2;
    	$seqfile = "seq_v" . $number;
    	@subselectoptions = ("-range", $seqencesdefinition[$count], $seqencesdefinition[$count+1], $inputfile, ">", $seqfile);

	# Modification DC - 05/05/2017
    	#makesystemcall("./vsubseqselect " . join(' ',@subselectoptions)); #DC
    	makesystemcall("vsubseqselect " . join(' ',@subselectoptions)); #LK - DC replaced vsubseqselect by vsubseqselect2
	#my ($hour,$min,$sec) = Now();
	if($logOption){
   		print LOG " vsubseqselect @subselectoptions\n"; # DC replaced vsubseqselect by vsubseqselect2 , [$hour:$min:$sec]
	}
   }
   return @seqencesdefinition;
}

#------------------------------------------------------------------------------
# transform output file to a structure
sub trans_data
{
  my $file = $_[0];
  my $elem_nbr = 0;
  my @repetitions = ();
  my @temp;
  my $line;
 open(FD, $file)  or die "Error in opening the input file : $!";
 while( <FD> ) 
 {
     next if /^(\s)*$/;  # skip blank lines
     chomp;              # remove trailing newline characters
		$line = $_;
		$line =~ s/>/> /g;
      @temp = split(/ +/, $line);
      if($temp[0] =~ /^\>/) 
      {
			$repetitions[$elem_nbr] = Rep->new();
			$repetitions[$elem_nbr]-> Length($temp[1]);
			$repetitions[$elem_nbr]-> Pos1($temp[2]);
     	} 
     else
     {
			$repetitions[$elem_nbr]-> DRseq($temp[0]);
			$elem_nbr++;
     }
  }
close(FD);
#my ($hour,$min,$sec) = Now();

if($logOption){
  print LOG "\n Getting results from vmatch and transform vmatch output file...\n"; #[$hour:$min:$sec]
}
  
return @repetitions;
}

#------------------------------------------------------------------------------
# find clusters and check DRs ---- 
sub write_clusters{

  my ($RefSeq,@rep) = @_;
  my @tabseq = find_clusters(@rep);

  #print "@tabseq\n"; #DC

  #my ($hour,$min,$sec) = Now();
  
  if($logOption){
  	print LOG "\n Find CRISPRs candidates and check DRs...\n"; #[$hour:$min:$sec]
  }

  my($count, $seqbeg, $seqend,$i,$j, $DR, $DRocc, $nbrcris,$DRlength);
  $i=0;
  $j=0;
  $nbrcris = 0;
  my $OneSpacerCris_nbr = 0;  # represents all 'hypothetical' CRISPRs, not only those having one spacer?
  my $modf=-1;
  my @modifcode =();
  # extract all the candidate sequences and modify the seq delimiters
  @tabseq = extractsequence(@tabseq);	

  # --- analyze the candidate DR in each extracted sequences
  for($count=1; $count<= $#tabseq/2+1; $count++)
  {
		my @DR_cand = ();	# store all the DR candidates
		my @DR_cand_occ = ();	# store the possible DRs after comparing the number of occurrences
		$seqbeg = $tabseq[2*$count-2]; #start position of the locus file
		$seqend = $tabseq[2*$count-1]; #end position of the locus file

		# -- store the length of the first DR candidat
		# --------------------
		$DR = $rep[$i]->DRseq;
		push(@DR_cand, $DR);

		my $occ_DR_max = DR_occ_rev($DR, @rep);
		
		push(@DR_cand_occ,$DR);
		$i++;
		# --------------------

		while(($i<= $#rep) && ($rep[$i]->Pos1 <= $seqend) )
		{
			$DR = $rep[$i]->DRseq;
			if(!check_DR($DR, @DR_cand))
			{
				push(@DR_cand, $DR);
				$DRocc = DR_occ_rev($DR, @rep);  #number of occurrences (rev comp included)
				# find the maximal number of occurrences (if equal, take the minimal length)
				if($occ_DR_max < $DRocc)
				{
					$occ_DR_max = $DRocc; 
					@DR_cand_occ=(); 
					push(@DR_cand_occ,$DR);
				}
				else
				{
			  		if($occ_DR_max == $DRocc)
			  		{
						if( $DR ne $DR_cand_occ[0] ){push(@DR_cand_occ,$DR);}

			  		}
				}
			}
			$i++;
		}

   	my $indexname = "seq_v".$count;
   	my $crisprfile = "crispr_result_".$count;   
	my $actual_path = getcwd(); # Addition DC - 07/2017
   	my $bestscore = 100000; 
   	my $bestindex = 0;
   	my $c=0;  # counter of the DRs consensus
   	my $TotDRerr=0;my $TotDRerr_best=0;

	## DC - 09/05/2017
	## Current repository AND Modifications chdir
	#my $currentRepositoryDC = getcwd(); #DC
	#print "CURRENT REP: $currentRepositoryDC\n"; #DC
	#chdir(".."); #DC
	#chdir(".."); #DC


   	while($c <= $#DR_cand_occ)
   	{
		my $score;
		# case we have more than one consensus for the DR
		if($#DR_cand_occ > 0)
		{
			# we have c DRfile and c crispr_result files
			$crisprfile = "crispr_result_".$count."_".$c;
		}	
		$DRlength = length($DR_cand_occ[$c]);
   		my $err = 0;  

		my @fuzznucoptions = ();
		
		if($betterDetectTruncatedDR){
			# DC better checking truncatedDR - begin
			my @charSeqDR = split(//, $DR_cand_occ[$c]);
			my $newDRpattern = "";
			for(my $k = 0; $k < $DRlength; $k++){
				if($k < (int($DRlength/2)) )
				{
					$newDRpattern .= $charSeqDR[$k];
				}
				else{
					$newDRpattern .= "n";
				}
			}
			#print "\nDR $c : $newDRpattern\n";
			$err = $betterDetectTruncatedDR; #int($DRlength/10);
			@fuzznucoptions = ("-sequence", $indexname, "-pattern", $newDRpattern, "-pmismatch", $err, "-outfile", $crisprfile);
			# DC better checking truncatedDR - end
		}
		else{
			$err = int($DRlength/$DRtrunMism); # DC - replaced "int($DRlength/3)" by "int($DRlength/3)"
   			@fuzznucoptions = ("-sequence", $indexname, "-pattern", $DR_cand_occ[$c], "-pmismatch", $err, "-outfile", $crisprfile);
		}
		push @fuzznucoptions, "-auto";
		
		makesystemcall("fuzznuc " . join(' ',@fuzznucoptions)); #DC
		#my $testCrisprFile = "Test_".$crisprfile;
		#makesystemcall("cp $crisprfile $testCrisprFile ");
		#my ($hour,$min,$sec) = Now();

		if($logOption){
			print LOG " fuzznuc " . join(' ',@fuzznucoptions)." \n"; #[$hour:$min:$sec]
		}

		if($#DR_cand_occ > 0)
		{
			#compare the crispr analysis based on the possible DRs
			($score,$TotDRerr) = compute_Crispr_score($crisprfile,length($DR_cand_occ[$c]));
			if($score <= $bestscore){$bestscore = $score; $bestindex = $c;$TotDRerr_best=$TotDRerr;}
			#print "Score = $score , Total DR err = $TotDRerr\n\n";
		}
		else
		{
			($bestscore,$TotDRerr_best) = compute_Crispr_score($crisprfile,length($DR_cand_occ[0]));
			#print "Score = $bestscore , Total DR err = $TotDRerr_best\n\n";
		}
		$c++;
   	} # boucle while
   	# take only the good DR seq and crispr_result seq (rename it to be good!)
   	if($#DR_cand_occ > 0)
   	{
   		$crisprfile = "crispr_result_".$count;
			rename($crisprfile."_$bestindex", $crisprfile);
			$DR = $DR_cand_occ[$bestindex];
   	}
   	else
	{
		$DR = $DR_cand_occ[0];
	}

   	my ($crisOK,$simDRs,$CrisprBeg, $CrisprEnd,$RefSpacersH,$nbspacers,$RefFalsSpacers) =
			Find_theCrispr($indexname,$DR,$count,$seqbeg,$seqend,$crisprfile);
	if($crisOK)
	{
		if($TotDRerr_best > $DRerrors){$crisOK = 0;} # DC - replace 0.2 by $DRerrors
	} 
   
   	# compute the number of crisprs in the whole sequence
   	$nbrcris = $nbrcris + $crisOK;
	
   	# if the current subsequence contains a crispr, it is stored
   	if($crisOK)
   	{
			my $Crispr_file = fill_in_crisprfile($simDRs,$RefSeq,$ResultDir, $nbrcris, $CrisprBeg, $CrisprEnd,$DR, $nbspacers, "spacers".$indexname,$RefFalsSpacers, %$RefSpacersH);

			#if($Crispr_file eq "ND") {$crisOK = 0;} # DC

			my @FalsSpacers = @$RefFalsSpacers;
			if($#FalsSpacers != -1){$modf=$#FalsSpacers;}
			if($nbspacers <= 1){$OneSpacerCris_nbr ++;} # DC replaced '$nbspacers <= 1 || $simDRs == 0' by '$nbspacers <= 3' or '$nbspacers <= 2'
			
			my $actual_path_after_fill_in = getcwd(); #DC
  			#print "ACTUAL PATH after fill_in_crisprfile !!!!: $actual_path_after_fill_in\n"; #DC
			
			#($hour,$min,$sec) = Now(); # DC - 07/2017
			
			if($logOption){
				print LOG "\n Nb spacers = $nbspacers , similarity DRs = $simDRs , Hypotheticals = $OneSpacerCris_nbr .\n"; #[$hour:$min:$sec]
			}


			if($modf != -1)
			{
				#DC - 09/05/2017
				#print "CHDIR located in Write Clusters : Path = $ResultDir / $RefSeq \n"; # DC
				# DC - 07/2017 - DANGEROUS MODIFICATION
				#($hour,$min,$sec) = Now(); # DC - 07/2017
				#$actual_path = getcwd();
				if($logOption){
					#print LOG "\n[$hour:$min:$sec] Actual path directory before CRISPR modification: $actual_path ...\n"; # DC - 07/2017
					#print "\n[$hour:$min:$sec] Actual path directory before CRISPR modification: $actual_path ...\n";
				}

		 		#chdir($ResultDir."/".$RefSeq); #
				chdir($ResultDir."/".$RefSeq); # Replace chdir($ResultDir."/".$RefSeq."/") by chdir($RefSeq."/")

			# DC - 07/2017 - End of DANGEROUS MODIFICATION

	    		my $hyp_cris_nbr = $OneSpacerCris_nbr;
		 		my $crisprs_nbr = $nbrcris;
		 		#--------------------------------------
		 		while($modf >=0)
		 		{
					#($hour,$min,$sec) = Now();
  
					if($logOption){
						print LOG "\n Modification would be performed (Nb good CRISPRs = $crisprs_nbr; Nb hypothetical = $hyp_cris_nbr)...\n"; #[$hour:$min:$sec]
					}

	  		 		if($crisprs_nbr - $hyp_cris_nbr >= 1)
	  		 		{
			  			my $s = 1;
			 		 	while($s<=$crisprs_nbr)
			 		 	{
							#DC - 05/2017 - replace $inputfile by $inputfileTmp (in the whole 'while' loop)
				     			#$inputfile = "$RefSeq"."_Crispr_".$s; # DC
							my $inputfileTmp = "$RefSeq"."_Crispr_".$s; # DC
				     			if(-e $inputfileTmp)
				     			{
								($s,$crisprs_nbr,$OneSpacerCris_nbr) =
									modify_files($inputfileTmp,$ResultDir,
										$RefSeq,$s,$crisprs_nbr,
										$OneSpacerCris_nbr,$modf);
									$nbrcris = $crisprs_nbr;
				     			}
				     			$s++;
			 		 	}
			 		}
	 	   		$modf--;
				}
			  #chdir($actual_path_after_fill_in); # DC - 07/2017 - IMPORTANT addition
			  #chdir($actual_path);
			  #if($logOption){
			  #	print LOG "\n[$hour:$min:$sec] Actual path directory after CRISPR modification: $actual_path ...\n"; # DC - 07/2017
			  #	print "\n[$hour:$min:$sec] Actual path directory after CRISPR modification: $actual_path ...\n";
			  #}

			}
			$modf = -1;
		
			if($logOption){
				print LOG "\n Actual path directory after CRISPR modification: $actual_path ...\n"; # DC - 07/2017 , [$hour:$min:$sec]
				#print "\n[$hour:$min:$sec] Actual path directory after CRISPR modification: $actual_path ...\n";
			}

			chdir($actual_path_after_fill_in); # DC - 07/2017 - IMPORTANT addition DC removed chdir($actual_path);
			

   	}
 } #boucle for

return ($nbrcris, $OneSpacerCris_nbr);
}

#------------------------------------------------------------------------------
# Function allowing ro modify files when a spacer is too long...
sub modify_files
{
  my($inputfile,$ResultDir,$RefSeq,$rank,$crisprs_nbr,$HYPCris_nbr,$modf)=@_;
  my(@lines,@temp,$nb_spacers,@nb_div, $newCrisNbr, @divi);

  #my ($hour,$min,$sec) = Now();
  my $actual_pathDC = getcwd(); #DC

  if($logOption){
  	print LOG " Modify CRISPRs files generated ($inputfile)...........\n Actual path: $actual_pathDC\n\n"; #[$hour:$min:$sec]
  }

  # Modification DC - 05/05/2017
  #my $dir = $ResultDir."/".$RefSeq."/"; #DC - 07/2017
  
  
  #print "ACTUAL PATH: $actual_pathDC\n"; #DC 
  my $dir = $ResultDir."/".$RefSeq."/"; # 05/2017  addition of "$ResultDir."/"."

  chdir($dir);

  open(FD, $inputfile)  or die "Error in opening the file:$inputfile!";
  @lines=<FD>;
  close(FD);

  $newCrisNbr = $crisprs_nbr;
  @temp = split(/:/, $lines[15]);
  my @temp2 = split(/ +/, $temp[1]);
  my $begPos = $temp2[1];
  my $endPos = $temp[2];
  @temp = split(/:/, $lines[16]);
  @temp2 = split(/ +/, $temp[2]);
  my $DRlength = $temp2[1];
  $nb_spacers = $temp[3];
  my $j = 20+2*$nb_spacers;

  if($lines[$j] =~ '#'){}
  else{
	# case we have divisions
	@temp = split(/:/, $lines[$j]);
	@nb_div = split(/ +/, $temp[1]);
	if($modf != -1) { @divi = @temp; shift @divi;shift @divi;pop @divi; }
		$newCrisNbr += scalar(@nb_div);
	# hypothese : there is only one division per CRISPR
	my $line = find_Sp_lines($j,$inputfile,$nb_div[0]);
	my($file_new, $file_old,$l1,$l2,$id,$CrisprBeg,$CrisprEnd,$nbspacers,$Spfile_old,$Spfile_new);
	$file_old = $inputfile;
	$Spfile_old = "Spacers_".$rank;

my ($nbsp1,$nbsp2);
# ----------file1-----------------------
$l1 = 20;
$l2 = $line ;
#print "rank = $rank inputfile : $inputfile\n-------------\n";
$id = $rank;
$Spfile_new = "Spacers_test_".$id;
$CrisprBeg = $begPos;
$CrisprEnd = $nb_div[0];
$CrisprEnd = $CrisprEnd."\n";
$nbsp1 = ($l2 - $l1)/2;
if($nbsp1>=1){
	#if($nbsp1 >= 2){
		$file_new = "$inputfile"."_test_".$id;
	#}    # DC replaced ">=" by ">"..... then ">= 3" replaced by ">= 2"
	#else{
	#	$file_new = "$inputfile"."_test_possible_".$id;$HYPCris_nbr++;
	#}
	
	$crisprs_nbr=create_file($file_new, $file_old,$l1,$l2,$id,$CrisprBeg,$CrisprEnd,$nbsp1,$modf,1,$crisprs_nbr,@divi);

	# DC - 05/05/2017
	#print " 1rst Calling create_spFile, DIR = $dir , ResultDir = $ResultDir\n";
	if($logOption){
  		print LOG " arguments of create_spFile : $dir,$ResultDir,$RefSeq,$Spfile_new, $Spfile_old, 0, $nbsp1\n"; #[$hour:$min:$sec]	
	}

	create_spFile($dir,$ResultDir,$RefSeq,$Spfile_new, $Spfile_old,0, $nbsp1);
}
else{$crisprs_nbr--;}

$l1 = $line + 2;
$l2 = $j ;
$id = $rank + 1;

$Spfile_new = "Spacers_test_".$id;

# get the wrong spacer length
my $sp = getSp($file_old,$l1);

$CrisprBeg = $sp - $DRlength;
$CrisprEnd = $endPos;
my $sp_prev = $nbsp1+1;
$nbsp2 = ($l2 - $l1)/2;
if($nbsp2>=1){
	#if($nbsp2 >= 2){ # DC replaced ">=" by ">"..... then ">= 3" replaced by ">= 2"
		$file_new = "$inputfile"."_test_".$id;
	#}   
	#else{
	#	$file_new = "$inputfile"."_test_possible_".$id;$HYPCris_nbr++;
	#}

	$crisprs_nbr =create_file($file_new, $file_old,$l1,$l2,$id,$CrisprBeg,$CrisprEnd,$nbsp2,$modf,2,$crisprs_nbr,@divi);

	# DC - 05/05/2017
	#print " 2nd Calling create_spFile, DIR = $dir , ResultDir = $ResultDir\n";
	create_spFile($dir,$ResultDir,$RefSeq,$Spfile_new, $Spfile_old, $sp_prev,$nbsp2);
}
else{
	$crisprs_nbr--;
}

$rank = $id;

# renaming the files #DC removed all "$dir/"
my $l;

if($nbsp1>=1 && $nbsp2>=1){
	for(my $k=$crisprs_nbr-1; $k>=$rank; $k--){
   		$l = $k + 1;
   		if(-e "$RefSeq"."_Crispr_".$k){
			rename("$RefSeq"."_Crispr_".$k, "$RefSeq"."_Crispr_".$l);
			rename("Spacers_".$k, "Spacers_".$l);
   		}
		#else{
		#	rename("$RefSeq"."_PossibleCrispr_".$k, "$RefSeq"."_PossibleCrispr_".$l);
		#}
	}
}
else{

	for(my $k=$crisprs_nbr-1; $k>=$rank; $k--){
   		$l = $k - 1;
   		if(-e "$RefSeq"."_Crispr_".$k){
			rename("$RefSeq"."_Crispr_".$k, "$RefSeq"."_Crispr_".$l);
			rename("Spacers_".$k, "Spacers_".$l);
   		}
		#else{
		#	rename("$RefSeq"."_PossibleCrispr_".$k, "$RefSeq"."_PossibleCrispr_".$l);
		#}
	}
}

$l = $rank - 1;

if($nbsp2>=1){
	if($nbsp1>=1){$l=$rank;}
	unlink("$RefSeq"."_Crispr_".$rank);
	if(-e "$inputfile"."_test_".$rank){
		rename("$inputfile"."_test_".$rank,"$RefSeq"."_Crispr_".$l);
	}
	#else{
	#	rename("$inputfile"."_test_possible_".$rank,"$RefSeq"."_PossibleCrispr_".$l);
	#}

	# ___spacers files ----------
 	$Spfile_new = "Spacers_test_".$rank; 
 	my $Spfile = "Spacers_".$l; 
 	rename("$Spfile_new", "$Spfile");
}

$l = $rank - 1;

if($nbsp1 >=1){ 
	unlink("$RefSeq"."_Crispr_".$l);
	if(-e "$inputfile"."_test_".$l){
		rename("$inputfile"."_test_".$l,"$RefSeq"."_Crispr_".$l);
	}
	#else{
	#	rename("$inputfile"."_test_possible_".$l,"$RefSeq"."_PossibleCrispr_".$l);
	#}

	# ___spacers files ----------

	$Spfile_new = "Spacers_test_".$l; 
	my  $Spfile = "Spacers_".$l; 
	rename("$Spfile_new", "$Spfile");
}

}

  return ($rank,$crisprs_nbr,$HYPCris_nbr);
}

#------------------------------------------------------------------------------
sub create_file{
 my($file_new, $file_old,$l1,$l2,$id,$CrisprBeg,$CrisprEnd,$nbspacers,$modf,$ordre,$crisprs_nbr,@divi) = @_;
 my $File_Content='';
 my @temp;

 #my ($hour,$min,$sec) = Now();

 if($logOption){
   print LOG "Create first CRISPRs files ...........\n"; #[$hour:$min:$sec] 
 }

  open(FD, $file_old)  or die "Error in opening the file:$file_old!";
   while(<FD>){
      if($. <= 14){$File_Content .= $_;}

      if($. == 15){
	$File_Content .= "# Crispr Rank in the sequence: $id\n";
        $File_Content .= "# Crispr_begin_position: $CrisprBeg\t Crispr_end_position: $CrisprEnd";
      }

      if($. == 17){
	@temp = split(/:/,$_);
	$File_Content .= $temp[0];
	$File_Content .= ":";
	$File_Content .= $temp[1];
	$File_Content .= ":";
	$File_Content .= $temp[2];
	$File_Content .= ":";
        $File_Content .= " $nbspacers\n";
        $File_Content .= "#=========================================================================\n";
        $File_Content .= "Spacer_begin_position\t Spacer_length\t Spacer_sequence\n";
      }
      if ( ($.>=$l1) && ($. < $l2) ){$File_Content .= $_;}
  }
  $File_Content .= "#=========================================================================\n";


 if(($#divi == -1)||($ordre == 1)){$File_Content .="########################################\n";}

else{
    #$File_Content .= "#Potentially false CRISPR; Spacers divisions:"; # DC
    #foreach(@divi){$File_Content .= $_;$File_Content .=":";}
    $File_Content .= "Spacers divisions:";
    foreach(@divi){$File_Content .= $_;$File_Content .=":";}
    $File_Content .="\n########################################\n";

  }

if($ordre == 2){$crisprs_nbr++;}


  close(FD);
  open WRITER,"> $file_new" or die "The file $file_new cannot be edited !\n";
  print WRITER $File_Content;

  if($logOption){
  	print LOG " File created: $file_new...........\n"; #[$hour:$min:$sec]
  }
  close(WRITER);

return $crisprs_nbr;
}
#####################################
#------------------------------------------------------------------------------
# Function allowing to create new spacer files when an array is modified
sub create_spFile{
 my ($dir,$ResultDir,$RefSeq,$Spfile_new,$Spfile_old,$sp_prev,$nbsp) = @_;
 $sp_prev = $sp_prev *2 +1;  # DC replaced 2 by 3 ?
 $nbsp = $sp_prev + $nbsp *2 -1; # DC replaced 2 by 3 ?

 #my ($hour,$min,$sec) = Now();
 # DC - 05/05/2017
 if($logOption){
 	print LOG "Old spacer file: $Spfile_old \n"; #DC , [$hour:$min:$sec] 
	print LOG " New spacer file: $Spfile_new \n"; #DC
	print LOG " Previous number of spacers: $sp_prev \n"; #DC
	print LOG " Number of spacers: $nbsp \n"; #DC

	#print "[$hour:$min:$sec] Old spacer file: $Spfile_old \n"; #DC
	#print "[$hour:$min:$sec] New spacer file: $Spfile_new \n"; #DC
	#print "[$hour:$min:$sec] Previous number of spacers: $sp_prev \n"; #DC
	#print "[$hour:$min:$sec] Number of spacers: $nbsp \n"; #DC
 }

 my $File_Content = "";

 open(FD, "$Spfile_old")  or die "Error in opening the file:$Spfile_old!"; # DC - 07/2017 - replace "$dir/$Spfile_old" by "$Spfile_old"
 
 # Modification DC - 05/05/2017
 #open(FD, "$Spfile_old")  or die "Error in opening the file:$Spfile_old!"; # DC replaced "$dir/$Spfile_old" by "$dir$Spfile_old", then by "$Spfile_old"


   while(<FD>){
	if( ($. >= $sp_prev) && ($. <= $nbsp) ){
	  $File_Content .= $_;
	}
   }
 close(FD);
 #Modification DC - 05/05/2017
 #open WRITER,"> $Spfile_new" or die "The file $Spfile_new cannot be edited !\n"; #DC
 open WRITER,"> $Spfile_new" or die "The file $Spfile_new cannot be edited !\n"; # ReDo - DC - 11/2017 removed "$dir/"

 print WRITER $File_Content;
 close(WRITER);
 
}
#------------------------------------------------------------------------------
sub getSp{
 #DC - 05/2017 - remove ','
 my($file_old,$l1,) = @_;
 #my($file_old,$l1,) = @_;
 my (@temp, $Sp);
  
  open(FD, $file_old)  or die "Error in opening the file:$file_old!";
   while(<FD>){
      if ($. == $l1){@temp = split(/ +/,$_);}
  }
  close(FD);
  $Sp = $temp[1];
#print "taille : $Sp\n";
return $Sp;
}
#------------------------------------------------------------------------------

sub find_Sp_lines{
 my($j,$inputfile,$nb_div) = @_;
 my $i=0;
 my $line;

 open(FD, $inputfile)  or die "Error in opening the file:$inputfile!";
 while(<FD>){
     if(($_ =~ $nb_div)&&($. != $j+1)){$line= $. ;} 
 }
 close(FD);
return $line;
}

#------------------------------------------------------------------------------
sub create_recap{
  #use Date::Calc qw(:all);
  my($RefSeq, $nbrcris, $OneSpacerCris_nbr, $ResultDir) = @_;
  my $File_Content;
  my $Crispr_report = $RefSeq."_CRISPRs"; 
  my $directory = $ResultDir;
  #my ($hour_recap,$min_recap,$sec_recap) = Now();

  if($logOption){
  	print LOG " Create first CRISPR(s) file(s) ...........\n"; #[$hour_recap:$min_recap:$sec_recap]
  }

  unless(-d $directory){ mkdir $directory or die "$0: I can not create the folder $directory: $!\n" }
  my $dir = $directory."/".$RefSeq;

  unless(-d $dir){ mkdir $dir or die "$0: I can not create the folder $dir: $!\n" }
  
  #my($year,$month,$day, $hour,$min,$sec) = Today_and_Now();
  my($year,$month,$day, $hour,$min,$sec) = (0,0,0,0,0,0);
  $File_Content  = "########################################\n";
  $File_Content .= "# Program: CRISPR Finder \n";
  $File_Content .= "# Author: Ibtissem GRISSA\n";
  $File_Content .= "# Rundate (GMT): $day/$month/$year $hour:$min:$sec\n";
  $File_Content .= "# Report file: $Crispr_report\n";
  $File_Content .= "########################################\n";
  $File_Content .= "#=======================================\n";
  $File_Content .= "# \n";
  $File_Content .= "# Sequence: $RefSeq \n";
  $File_Content .= "# Number of CRISPR-like structures : $nbrcris\n";
  $File_Content .= "# Number of questionable structures : $OneSpacerCris_nbr\n";
  $File_Content .= "# \n";
  $File_Content .= "#=======================================\n";
  chdir($dir);
  open WRITER,"> $Crispr_report" or die "The file $Crispr_report cannot be edited !\n";
  print WRITER $File_Content;
  close WRITER;
}
#------------------------------------------------------------------------------
# copy spacers file of crispr_$refseq_$i
# this file will be edited by user (button)
sub copy_spacers{

  my($SpacersFile,$dir,$id) = @_; 
  $dir = $dir."/Spacers_$id";

  #print "j'ai copiÈ $SpacersFile dans $dir";
  #my ($hour,$min,$sec) = Now();
  
  if($logOption){
  	print LOG " Copy Spacers files in $dir...........\n"; #[$hour:$min:$sec]
  }

  my $cop = "cp -R $SpacersFile $dir";

  ## DC- 09/05/2017
  #my $cop = "cp $SpacersFile $dir"; #DC
  
  if($logOption){
  	print LOG " $cop ...........\n"; #[$hour:$min:$sec]
  }

  system($cop);
}
#------------------------------------------------------------------------------
# fill in the file crispr_$refseq_$i
# this file will be used for the database storage
sub fill_in_crisprfile
{
  #use Date::Calc qw(:all);

  my(	$simDRs,$RefSeq,$ResultDir,$id,$CrisprBeg,$CrisprEnd,
  	$DRseq,$nbspacers,$SpacersFile,$RefFalsSpacers,%spacersPos) = @_; 

  #my ($hour_crispr,$min_crispr,$sec_crispr) = Now();

  if($logOption){
  	print LOG " Create first CRISPRs files ...........\n"; #[$hour_crispr:$min_crispr:$sec_crispr]
  }

  my ($Crispr_file,$File_Content,$Bpos,$i);
  my $DRlength = length($DRseq);
  #my($year,$month,$day, $hour,$min,$sec) = Today_and_Now();
  my $directory = $ResultDir;
  if(-e "$directory"){}else{mkdir($directory);}
  my $dir = $directory."/".$RefSeq;
  if(-e "$dir"){}else{mkdir($dir);}

  #if(($nbspacers > 1) && ($simDRs ==1)) # DC replaced '$nbspacers >= 2' by '$nbspacers > 3' AND '(($nbspacers > 3) && ($simDRs ==1))' replaced by ($nbspacers > 3) 

  #if($nbspacers >= 1)
  #{ # the replaced '> 3' by '>=3' then by '>=2' then by '>2'
  	$Crispr_file = "$dir/$RefSeq"."_Crispr_".$id;
  #}
  #else
  #{
  #	$Crispr_file = "$dir/$RefSeq"."_PossibleCrispr_".$id;
  #}
  $nbspacers ++;
  $File_Content  = "########################################\n";
  $File_Content .= "# Program: CRISPR Finder\n";
  $File_Content .= "# Author: Ibtissem GRISSA\n";
  $File_Content .= "# Rundate (GMT): ND\n";  #$day/$month/$year $hour:$min:$sec\n";
  $File_Content .= "# Report_file: $Crispr_file\n";
  $File_Content .= "########################################\n";
  $File_Content .= "#=======================================\n";
  $File_Content .= "# \n";
  $File_Content .= "# Sequence: $RefSeq \n";
  $File_Content .= "# Description: ".$seq->desc()."\n";
  $File_Content .= "# Length: ".$seq->length()."\n";
  $File_Content .= "# Id: ".$seq->id()."\n";
  $File_Content .= "#\n";
  $File_Content .= "#=========================================================================\n";
  $File_Content .= "# Crispr Rank in the sequence: $id\n";
  $File_Content .= "# Crispr_begin_position: $CrisprBeg\t Crispr_end_position: $CrisprEnd\n";
  $File_Content .= "# DR: $DRseq\t DR_length: $DRlength\t Number_of_spacers: $nbspacers\n";
  $File_Content .= "#=========================================================================\n";
  $File_Content .= "Spacer_begin_position\t Spacer_length\t Spacer_sequence\n";

  copy_spacers($SpacersFile, $dir,$id);
  open(FD, $SpacersFile)  or die "Error in opening the input file : $SpacersFile!";
  my @lines=<FD>;
  close(FD);

  # DC - 05/2017 - Detect identical sapcers
  my %hashCountSpacers = ();  # DC - a hash table to count spacers
  my $duplicatedSpacers = 0; # DC - boolean variable allowing to detect duplicated spacers
  my @tempChainSpacer = (); # DC - table containing duplicated spacers

  $i=1;
  foreach my $k (reverse (sort { $b <=>$a } keys %spacersPos) )
  {
	$Bpos = " "x(20 - length($k));
	$File_Content .= "$Bpos $k";
	$Bpos = ' 'x(12 - length($spacersPos{$k}));
	$File_Content .= "\t $Bpos $spacersPos{$k}";
	$File_Content .= "\t $lines[$i]\n";

	$hashCountSpacers{$lines[$i]}++; # DC - increment the hash table to calculate occurrence of each spacer

	$i +=2;
  }
  $File_Content .= "#=========================================================================\n";

  # DC - 05/2017 - Show repeated spacers
  foreach my $k2 (keys(%hashCountSpacers)) {
	if($hashCountSpacers{$k2} > 1){
		$duplicatedSpacers = 1;
		push(@tempChainSpacer,$k2);
  		
	}
  } # DC
  
  #if($duplicatedSpacers == 1){
  #	$File_Content .= "#Potentially false CRISPR; Repeated spacers (n=".($#tempChainSpacer+1).").";
  #  	$File_Content .= "\n"; 
  #} # DC

  
  my @FalsSpacers = @$RefFalsSpacers;
  if($#FalsSpacers == -1){$File_Content .="########################################\n";}
  else
  {
    #$File_Content .= "#Potentially false CRISPR; Long Spacers :"; # DC
    #foreach(@FalsSpacers){$File_Content .= $_;$File_Content .=":";}
    $File_Content .= "Spacers divisions:";
    foreach(@FalsSpacers){$File_Content .= $_;$File_Content .=":";}

    $File_Content .="\n";
  }
  $File_Content .="########################################\n";

  open WRITER,"> $Crispr_file" or die "The file $Crispr_file cannot be edited !\n";
  print WRITER $File_Content;
  close WRITER;

  if($logOption){
  	print LOG "[ Create file: $Crispr_file ...........\n"; # $hour_crispr:$min_crispr:$sec_crispr]
  }

  return $Crispr_file;
}

#------------------------------------------------------------------------------
# compute an error score for a DR matches in the crispr file
# used in the comparison of two DRs
sub compute_Crispr_score
{
  my($crisprfile,$DRlength) = @_;
  my $score = 0;
  my(@lines, @temp, $penal, $TotERR,$nb);
  my $troncated = 0;  
  $TotERR =0; $nb=0;
 
  #my %testHash = (); # DC test
  #my $fileTest = "/home/davidcouvin/Applications/CRISPRCasFinder/HASH_repeats.txt"; # DC test

#print "<br>";
  open(FD, $crisprfile)  or die "Error in opening the input file : $crisprfile!";
  @lines=<FD>;
  close(FD);
  my $count = 19;
  if($#lines <= 5){$score = 100000;}	#added in 15/11/06	
  for($count=12; $count < scalar(@lines)-1; $count++)
  {
     $nb++;
     next if $lines[$count] =~ /^#/;
     next if $lines[$count] =~ /^\s+$/;
     $lines[$count] = trim($lines[$count]);
     next if $lines[$count] =~ /^Start/;
     @temp = split(/ +/, $lines[$count]);
     #if($temp[0]=~""){shift @temp;}

     #DC
     #open(HASH, ">>$fileTest")  or die "Error in opening the input file : $fileTest!";
     #print HASH "\nSEQ in function = $temp[$#temp]\n"; 
     #close(HASH);

     if($force){
	
	#print "\n******\nMismatch = $temp[$#temp-1] ; SEQ in function = $temp[$#temp]\n";
	if($temp[$#temp] =~ /^G/){ $score = $score -2; } #DC code to modify score (when e.g. G is encountered at the beginning) 
	if(length($temp[$#temp]) >= 30){ $score = $score -2; } #DC check DR size
	if($temp[$#temp] =~ /AA.$/){ $score = $score -2; } #DC code to modify score (when e.g. AAN is encountered at the end)
     }

     if($temp[$#temp-1] ne "."){
	$TotERR = $TotERR + $temp[$#temp-1];
	$penal = 1 + $temp[$#temp-1]/$DRlength;
        #print "$count penal: $penal <br>";
	$score = $score + $penal;
	if($count == 19){$troncated = $penal;}
	if($count == scalar(@lines)-4){
		if($penal > $troncated){
			$troncated = $penal;
		}
	}
     }
  }
  #print "tronc : $troncated <br>";
  $score = $score - $troncated;
  $TotERR = $TotERR/($nb*$DRlength);
  #print "<br>$crisprfile score : $score __ toterr : $TotERR <br>";
  #$ans = <STDIN>;
  return ($score,$TotERR);
}
#------------------------------------------------------------------------------
# study a given seq_v according to a specific candidate DR
# 26/10/06 use of clustalw to align spacers 
# March, 2018, ClustalW has been replaced by Muscle
sub Find_theCrispr
{
   my($indexname,$DR,$count,$seqbeg,$seqend,$crisprfile)= @_;
   my(@spacers, $goodcrispr,$CrisBeg,$CrisEnd,$refFalsSpacers,$simDRs);
   my $DRlength = length($DR);
   
   ($simDRs,$refFalsSpacers, @spacers) = definespacers($DRlength, $crisprfile,$indexname);
   my @FalsSpacers = @$refFalsSpacers;
   if($#FalsSpacers >=0){for(my $h =0; $h<= $#FalsSpacers; $h++){$FalsSpacers[$h] += $seqbeg;}}
   my %spacersH = trans_struct_hach($seqbeg,@spacers);
   
   #my ($hour,$min,$sec) = Now();
   
   if($logOption){
   	print LOG " Check short CRISPRs arrays and check spacers alignment using Muscle....\n"; #[$hour:$min:$sec]
   }

   $goodcrispr = 0;
   my $o = $#spacers;

   #print "\nDR: $DR\n";
   if($minNbSpacers <= 0){ $minNbSpacers = 1;}

   if( $#spacers >= ($minNbSpacers - 1) ) # MINNBSPACERS
   {
	if($#spacers == 0)
	{
		$goodcrispr = check_short_crispr($DR, "spacers".$indexname);
	}
	else
	{
		$goodcrispr = checkspacersAlignMuscle("spacers".$indexname);
	}
	if($goodcrispr)
	{
		$CrisBeg=0; $CrisEnd=0;
		$CrisBeg = $seqbeg + $spacers[0]->Pos1 - $DRlength;
		$CrisEnd = $seqbeg + $spacers[$#spacers]->Pos1 + $spacers[$#spacers]->Length +$DRlength -1;
	}
  }
  my $r=$#FalsSpacers;
  
  return($goodcrispr,$simDRs,$CrisBeg,$CrisEnd,\%spacersH,$#spacers,\@FalsSpacers);
}

#------------------------------------------------------------------------------
# transform the old @spacers structure to a hachage containing as keys the
# begin position of a spacer (according to the whole sequence) and as values
# the corresponding length
sub trans_struct_hach
{
   my($seqbeg,@structspacers) = @_;
   my %spacers = ();
   my ($i,$pos, $len);
   for($i=0; $i<= $#structspacers; $i++)
   {
	$pos = $structspacers[$i]-> Pos1 + $seqbeg;
	$len = $structspacers[$i]-> Length;
	$spacers{$pos} = $len;
   }
   return %spacers;
}

#------------------------------------------------------------------------------
# this function is used to check, in the case of a unique spacer, the identity 
# between the DR and the CRISPR 
# use of program needle of emboss to do the alignment

sub check_short_crispr
{
  my($DR, $spacerfile) = @_;
  my($needleoptions, $goodcrispr);

  # edit DR file
  open WRITER,"> DR" or die "The file cannot be edited !\n";
  print WRITER "> DR\n";
  print WRITER $DR;
  print WRITER "\n";
  close WRITER;

  # needle program
  #$needleoptions = "needle -auto -asequence DR -bsequence $spacerfile -gapopen 10.0 -gapextend 0.5 -outfile alignDR_Spacer.needle"; 
  $needleoptions = "needle -auto -asequence DR -bsequence $spacerfile -gapopen 10.0 -gapextend 0.5 -outfile alignDR_Spacer.needle"; # DC
  makesystemcall($needleoptions);

  #my ($hour,$min,$sec) = Now();
  
  if($logOption){
  	print LOG " $needleoptions\n"; #[$hour:$min:$sec]
  }

  $goodcrispr = similarity_needle();
#print "goodneedle = $goodcrispr <br>";
  return $goodcrispr;
}
#------------------------------------------------------------------------------
# analyze needle file to conclude if repetitivity is significant
sub similarity_needle{

  open(FD, "alignDR_Spacer.needle")  or die "Error in opening the input file : alignDR_Spacer.needle!";

  #my ($hour,$min,$sec) = Now();
  # find the gap value
  while(<FD>){
	if($_ =~ "Gaps:" ){
	  my @temp = split(/ +/, $_);
	  #print "SEEEEEE THIS LINE: $_ \n";
	  my $score = $temp[3];
	  $score =~ s/\(//;
	  $score =~ s/%\)//;
	
          if($logOption){
  		print LOG " Running function similarity_needle... score is: $score\n"; #[$hour:$min:$sec]
  	  }	
          if(!$score){ $score=0; }  

	  if($score > 50.0){return 1;}else{return 0;}
	}
  }
  close(FD);
}

#------------------------------------------------------------------------------
# checks if the DR is already a candidate 
sub check_DR{
  my ($DR, @DR_cand) = @_;
  my $i=0;
  my $stop = 0;
  
  while(($i <= $#DR_cand) && ($stop==0) ){
	if($DR eq $DR_cand[$i]){$stop = 1;};
	$i++;
  }
return $stop;
}
#------------------------------------------------------------------------------
# counts the number of occurrences of the input DR in the analyzed genome (the reverse complement is counted also)
sub DR_occ_rev
{
  my ($DR, @rep) = @_;
  my $DRocc = 0;
  my $count=0;
  my($seqDR, $revDR);

  $seqDR = Bio::Seq->new(-seq => $DR);
  $revDR = $seqDR->revcom();
  for($count=0; $count<=$#rep; $count++){
	if(($DR eq $rep[$count]->DRseq) || ($revDR->seq()eq$rep[$count]->DRseq) ){
		$DRocc++;
	}
  }
  return $DRocc;
}

#------------------------------------------------------------------------------
#find clusters and store their start and end positions
# @tabseq = find_clusters(@rep);
sub find_clusters		# find all possible clusters (even the shortest)
{
 my @repetitions = @_;
 my @tabseq = ();
 my $count = 0;
 my $nb_clust = 0;	# nb_clust = 2*nbr_clusters

 #my ($hour,$min,$sec) = Now();
 
 if($logOption){
 	print LOG " Find clusters and store their start and end positions\n"; #[$hour:$min:$sec]
 } 

#print "nbr de rep = $#repetitions --\n";
$tabseq[$nb_clust] = $repetitions[$count]->Pos1;	# start position of the first cluster
  while($count < $#repetitions)
  {
    if(!compare_clusters($repetitions[$count]->Pos1,$repetitions[$count+1]->Pos1))
    {
     # a new cluster is found
     if($repetitions[$count]->Pos1 !=$tabseq[$nb_clust])
     {
	$nb_clust++; 
	$tabseq[$nb_clust] = $repetitions[$count]->Pos1;
      }
      else
      {
	$nb_clust++; 
	$tabseq[$nb_clust] = $repetitions[$count]->Pos1 + 60;

      }

     $nb_clust++; 
     $tabseq[$nb_clust] = $repetitions[$count+1]->Pos1;
	
     }
     $count++;
   }


   if($#tabseq % 2 == 0)
   {
		if( $repetitions[$count]->Pos1 !=$tabseq[$nb_clust] )
		{
			$nb_clust++; 
			$tabseq[$nb_clust] = $repetitions[$count]->Pos1;
		}
		else
		{
			$nb_clust++; 
			$tabseq[$nb_clust] = $repetitions[$count]->Pos1 + 60;
		}
   }
   return @tabseq;
}

#------------------------------------------------------------------------------
sub compare_clusters
{
  my($element1,$element2) = @_;
  if( ($element2 >= $element1 - 1500) && ($element2 <= $element1+1500) )
  {
    return 1;
  }else
  {
    return 0;
  }
}

#------------------------------------------------------------------------------
sub checkdbfile
{
  my($inputfile,$prjfile) = @_;
  my $dbfile = '';

  if(-e $prjfile)
  {
    unless(open(PRJFILEPTR,$prjfile))
    {
      return 0;
    }
    while(my $line = <PRJFILEPTR>)
    {
      if($line =~ /^dbfile=(\S+) (\d+)/)
      {
        if($dbfile eq '')
        {
          $dbfile = $1;
          my $dbfilesize = $2;
          if($dbfile eq $inputfile)
          {
            if(my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,
                   $atime,$mtime,$ctime,$blksize,$blocks) = stat($dbfile))
            {
              if($size == $dbfilesize)
              {
                close PRJFILEPTR;
                return 1;
              }
            }
          }
        }
      }
      close PRJFILEPTR;
      return 0;
    }
  }
  return 0;
}
#------------------------------------------------------------------------------
sub makesystemcall
{
  my ($argstring) = @_;
  my @args = split(' ',$argstring);
  my($retcode) = system($argstring);
  $retcode = $? >> 8;
  if($retcode ne 0)
  {
    print STDERR "failure: \"$argstring\", errorcode $?\n";
    exit 1;
  }
  #print STDERR "# $argstring\n"; #skip verbose mode
}

sub callmkvtree
{
  my ($inputfile,$indexname) = @_;
  if(not (checkdbfile($inputfile,"${indexname}.prj")))
  {
    makesystemcall("mkvtree -db $inputfile " .   # LK - DC replaced mkvtree by mkvtree2
                   "-dna -pl -lcp -suf -tis -ois -bwt -bck -sti1");
    #my ($hour,$min,$sec) = Now();

    if($logOption){
    	print LOG " mkvtree -db $inputfile -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1\n"; # DC replaced mkvtree by mkvtree2, [$hour:$min:$sec]
    }

  }
}

#-------------------------------------------------------------------------------
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
#------------------------------------------------------------------------------

sub printhelpall
{
print <<HEREDOC;
Name:
  $0 standalone version $version

Synopsis:
  A perl script to identify CRISPR arrays and associated cas genes in DNA sequences

Usage:
  perl $0 <filename.fasta>
  
  OR

  perl $0 [options] -in <filename.fasta>

  --Please note <filename.fasta> must be in Fasta format. Please also note that when several options are called, the option "-in or -i" must precede the input FASTA file. 

General:
  -help or -h		This help

  -version or -v	The current version of the program will be displayed

Other options:

  [Input/Output and -so]
  -in or -i [XXX]	Input Fasta file (with extensions: .fasta, .fna, .mfa, .fa, .txt)

  -outdir or -out [XXX]	Output directory (if users do not use this option, a default directory will be created with the date and time)

  -keepAll or -keep	Option allowing to keep secondary folders/files (Prodigal/Prokka, CasFinder, rawFASTA, Properties) (default: $keep)
 
  -LOG or -log		Option allowing to write LOG files (default: $logOption)

  -HTML or -html	Option allowing to display results as a static HTML web page (default value: $html). The web page created (index.html) will be dependent of a CSS file (crispr.css)

  -copyCSS [XXX]	Option allowing to copy provided CSS file into "Visualization" repository if option -HTML is set (default: '$cssFile')

  -soFile or -so [XXX]	Option allowing to use the shared object file if it is not present in current directory (default: '$so')
  
  -quiet or -q	Option allowing to run the program quieter (default: $quiet)

  -faster or -fast	Option allowing to run the program faster (default value: $fast)

  -minSeqSize or -mSS [XXX]	Option allowing to fix a sequence minimum size to search CRISPR and Cas systems (default: $seqMinSize)


  [Detection of CRISPR arrays]
  -mismDRs or -md [XXX]	Percentage mismatchs allowed between DRs (default: $DRerrors)

  -truncDR or -t [XXX]	Percentage mismatchs allowed for truncated DR (default: $DRtrunMism)

  -minDR or -mr [XXX]	Minimal size of DRs (default: $M1)

  -maxDR or -xr [XXX]	Maximal size of DRs (default: $M2)

  -minSP or -ms [XXX]	Minimal size of Spacers (default: $S1)

  -maxSP or -xs [XXX]	Maximal size of Spacers (default: $S2)

  -noMism or -n	Option used to do not allow mismatches (default value is $mismOne when this option is not called. i.e. mismatches are allowed by default)

  -percSPmin or -pm [XXX]	Minimal Spacers size in function of DR size (default: $Sp1)

  -percSPmax or -px [XXX]	Maximal Spacers size in function of DR size (default: $Sp2)

  -spSim or -s [XXX]	Maximal allowed percentage of similarity between Spacers (default: $SpSim)

  -DBcrispr or -dbc [XXX]	Option allowing to use a CSV file of all CRISPR candidates contained in CRISPRdb (from last update) (default: 'supplementary_files/CRISPR_crisprdb.csv')

  -repeats or -rpts [XXX]	Option allowing to use a consensus repeats list generated by CRISPRdb (default: 'supplementary_files/Repeat_List.csv')

  -DIRrepeat or -drpt [XXX]	Option allowing to use a file containing repeat IDs and orientation according to CRISPRDirection (default: 'supplementary_files/repeatDirection.tsv')

  -flank or -fl [XXX]	Option allowing to set size of flanking regions in base pairs (bp) for each analyzed CRISPR array (default: $flankingRegion)

  -levelMin or -lMin [XXX]	Option allowing to choose the minimum evidence-level corresponding to CRISPR arrays we want to display in Crisprs_REPORT file (default: $levelMin)

  -classifySmallArrays or -classifySmall or -cSA	Option allowing to change evidence level status of small arrays (evidence-level 1) having the same consensus repeat as an evidence-level 4 array (default value=$classifySmall)
  
  -forceDetection or -force	Option allowing to force/foster detection of specific CRISPR arrays (default: $force)

  -fosterDRLength or -fDRL [XXX]	Option allowing to foster a specific repeat's length when option '-force' is set (default: $fosteredDRLength)

  -fosterDRBegin or -fDRB [XXX]	Option allowing to foster a specific repeat's beginning when option '-force' is set (default: '$fosteredDRBegin'), regular expressions could be considered

  -fosterDREnd or -fDRE [XXX]	Option allowing to foster a specific repeat's ending when option '-force' is set (default: '$fosteredDREnd'), regular expressions could be considered 

  -MatchingRepeats or -Mrpts [XXX]	Option allowing to specify a query file containing repeats to be matched (default: '$repeatsQuery')

  -minNbSpacers or -mNS [XXX]	Option allowing to specify the minimum number of spacers required for each array (default: $minNbSpacers)

  -betterDetectTrunc or -bDT	Option allowing to better detect the truncated DR (default: $betterDetectTruncatedDR)

  -PercMismTrunc or -PMT [XXX]	Option allowing to set the percentage of allowed mismatches in the half of the truncated DR (default: $percentageMismatchesHalfDR)


  [Detection of Cas clusters]
  -cas or -cs	Search corresponding Cas genes using Prokka (default kingdom: "$kingdom") and MacSyFinder (default: $launchCasFinder)

  -ccvRep or -ccvr	Option used to write the CRISPR-Cas vicinity report (CRISPRs and Cas) if option -cas is set (default: $writeFullReport)

  -vicinity or -vi [XXX]	Option used to define number of nucleotides separating a CRISPR array from its neighboring Cas system (default: $vicinity)
  
  -cpuMacSyFinder or -cpuM [XXX]	Option allowing to set number of CPUs to use for MacSyFinder (default: $cpuMacSyFinder)

  -rcfowce	Option allowing to run Casfinder only when any CRISPR exists (default: $rcfowce) (set if -cas is set)

  -definition or -def [XXX]	Option allowing to specify CasFinder definition (if option -cas is set) to be more or less stringent (default: '$definition')

  -gffAnnot or -gff [XXX]	Option allowing user to provide an annotation GFF file (if options -cas and -faa are set) (default: '$userGFF') 

  -proteome or -faa [XXX]	Option allowing user to provide a proteome file '.faa' (if options -cas and -gff are set) (default: '$userFAA')

  -cluster or -ccc [XXX]	Option allowing to constitute clusters or groups of CRISPR or Cas systemes given a determined threshold e.g. 20kb (default: $clusteringThreshold) (set if -cas is set)

  -getSummaryCasfinder or -gscf	Option allowing to get summary file of Cas-finder (MacSyFinder) and copy it to TSV repository (default: $gscf)

  -geneticCode or -gcode [XXX]	Option allowing to modify the genetic code (translation table) for CDS annotation (default: $genCode)

  -metagenome or -meta	Option allowing to analyze metagenome (default: $metagenome)


  [Use Prokka instead of Prodigal (default option)] Prokka (https://github.com/tseemann/prokka) must be installed in order to use following options
  -useProkka or -prokka	Option allowing to use Prokka instead of Prodigal (default: $useProkka)

  -cpuProkka or -cpuP [XXX]	Option allowing to set number of CPUs to use for Prokka (default: $cpuProkka)

  -ArchaCas or -ac	same option as -cas using "Archaea" as default kingdom instead of "Bacteria" (default: $launchCasFinder). Option to be used when -prokka is used.


Options waiting for a given parameter (filename, text, or number) are followed by symbols "[XXX]". Other options could be considered as boolean (yes or no, 1 or 0).

#####################################################

The input file should meet these constraints:

- the file name must not contain any special characters or spaces, 
- the file name must not contain multiple dots (an acceptable file name is e.g. "multifasta.fna"),
- the sequence must be identified/named (the ID follows character ">", and a description could be added after a space character),
- the ID should not contain special characters such as "\$%" or multiple dots,
- The "|" sign can be used as field separator,
- the file must contain nucleotides (not amino acids),
- the file could contain several sequences in FASTA format,
- each ID must be unique,
- the ID and the file name must not be too long,
- the ID will be used for output.

Examples:
(1): perl $0 test.fasta
In this example, your result folder will be in the directory named: "Result_test"

(2): perl $0 -in test.fasta -md 20 -t 33.3 -mr 23 -xr 55 -ms 25 -xs 60 -pm 0.6 -px 2.5 -s 60

(3): perl $0 -in genomes100.fna -drpt supplementary_files/repeatDirection.tsv -rpts supplementary_files/Repeat_List.csv -cs -ccvr -dbc supplementary_files/CRISPR_crisprdb.csv -cf CasFinder-2.0.3 -html

(4): perl $0 -in sequence.fasta -cas -log -out RES_Sequence -ccc 20000 -ccvRep -keep -html -rcfowce -def S -cpuM 4 -copyCSS supplementary_files/crispr.css -cf CasFinder-2.0.3

(5): perl $0 -in sequence.fasta -cas -log -out RES_Sequence -force -so path/to/sel392v2.so

HEREDOC
}

sub printversion
{
    my $programname = shift @_;
    print "This is $programname, version $version,\n";
    print "a perl script to identify CRISPR arrays and associated cas genes from DNA sequences\n";
}

sub printhelpbasic
{
    my $programname = shift @_;
    print "Usage: $programname [options] [-in] <filename>\n";
    print "type -version or -v to see the current version\n";
    print "type -help or -h for help\n";
}

#!/usr/bin/env bash

vmatch="vmatch-2.3.0-Darwin_i386-64bit"
vmatch_url="http://vmatch.de/distributions/${vmatch}.tar.gz"

macsyfinder="macsyfinder-1.0.5"
macsyfinder_url="https://dl.bintray.com/gem-pasteur/MacSyFinder/${macsyfinder}.tar.gz"

casfinder="CasFinder-2.0.2"

prokka="prokka-1.12"
prokka_url="http://www.vicbioinformatics.com/${prokka}.tar.gz"

PREFIX='/usr/local/'
CURDIR=$(pwd)

stamp=$(date '+%Y-%m-%d-%H:%M:%S')
LOGFILE=${CURDIR}/installer.${stamp}.log
crispr_perllib="${PREFIX}libexec/crisprcas"

################### script below SHOULD NOT BE EDITED ########


####################
# helper functions #
####################

command_exists () {
    command -v "$1" > /dev/null 2>&1;
}

mac_install () {
    VOLUME=`hdiutil attach $1 | grep Volumes | awk '{print $3}'`
    cp -rf "${VOLUME}/*.app" "/Applications/"
    hdiutil detach ${VOLUME}
}

download () {
    PN="$1"
    url="$2"
    echo "Downloading ${url} ..."
    # -L follow redirection
    # -O to redirect output in file named as on server
    curl -L -O --silent "${url}"
    if [ $? -ne 0 ];
    then
        echo "ERROR: Download ${PN} from ${url} failed" | tee -a ${LOGFILE}
        echo "Download it manually and put the archive (as is) in this directory" | tee -a ${LOGFILE}
        exit 1
    fi
}

perl_pack_exists () {
    pack="$1"
    perllib="${crispr_perllib}/lib/perl5"
    PERL5LIB="${perllib}:${PERL5LIB}" perl -e "use ${pack} ;" > /dev/null 2>&1 ;
}

file_owner (){
    ls -ld $1 | awk '{print $3}'
}

file_group (){
    ls -ld $1 | awk '{print $4}'
}

brew_install () {
    pack="$1"
    cmd="$2"
    if ! command_exists ${cmd} ;
    then
        echo "Installation of ${pack}..." | tee -a ${LOGFILE}
        brew install ${pack} | tee -a ${LOGFILE}
        if [ $? -ne 0 ];
        then
            echo "ERROR: brew cannot install ${pack}, see ${LOGFILE} for details." | tee -a ${LOGFILE}
            exit 3
        fi
    fi
}



#################
# brew packages #
#################

# brew_packages is "like" an associative array
# (thanks to mac os to run very old bash which not support hash)
# keys are name of brew package
# value is the name of a command to test if it already exists

brew_packages=( "gcc:gcc"
                "wget:wget"
                "hmmer:hmmsearch"
                "blast:blastn"
                "cpanminus:cpanm"
                "emboss:seqret"
                "muscle:muscle"
                "parallel:parallel"
                "clustal-w:clustalw2"
                "prodigal:prodigal"
               )


#################
# perl packages #
#################
perl_packages="Try::Tiny Test::Most JSON::Parse Date::Calc Class::Struct Bio::DB::Fasta File::Copy Bio::Seq Bio::SeqIO"
perl_packages_forced="Bio::Tools::Run::Alignment::Clustalw Bio::Tools::Run::Alignment::Muscle"

####################
# Argument parsing #
####################

dependencies="Homebrew, java8, Bioperl, perl-tiny, perl-most, perl-json, perl-calc,\
 vmatch, macsyfinder , prokka, tbl2asn\n\
 ${brew_packages[@]%%:*} ${perl_packages} ${perl_packages_forced}"

echo -e "\nCRISPRCasFinder depends on: \n\
 ${dependencies} \n\
If they are not already installed this script will installed them."

while true; do
    read -p "Do you want to proceed (Y/N)?" yn
    case $yn in
        [Yy]* ) break ;;
        [Nn]* ) exit 0;;
        * ) echo "Please answer yes or no.";;
    esac
done


#########################
# Homebrew installation #
#########################
if ! command_exists brew ;
then
    echo "Homebrew installation ..." | tee -a ${LOGFILE}
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

    dir_2_check="bin etc include lib sbin share var Frameworks"
    my_groups=$(groups)
    for dir in ${dir_2_check};
    do
        path="/usr/local/${dir}"
        owner=$(file_owner ${path})
        group=$(file_group ${path})
        if [ ! -w "${path}" ];
        then
            echo "${path} is not writable brew need it"
            if [ ${USER} -eq ${owner} ]; then
                echo "Add permission to write for ${USER} on ${path}"
                sudo chmod u+w "${path}"
            elif [[ " ${my_groups[@]} " =~ " ${group}] " ]];
            then
                echo "Add permission to write for group \'${group}\' on ${path}"
                sudo chmod g+w "${path}"
            else
                echo "${path} is not writable."
                echo "brew need to have permission to write in /usr/local/${dir_2_check}."
                echo "Fix this issue and run this script again."
                exit 3
            fi
        fi
    done
else
   echo "Found installed Homebrew, skip installation."
fi

brew update

taps=$(brew tap)
tap="brewsci/science"
if [[ ! " ${tap[@]} " =~ " ${tap} " ]];
then
    echo "Adding tap ${tap}" | tee -a ${LOGFILE}
    brew tap ${tap}
    if [ $? -ne 0 ];
    then
        echo "ERROR: in adding ${tap} in homebrew" | tee -a ${LOGFILE}
        exit 3
    fi
else
    echo "${tap} .... found."
fi

brew update



#################
# brew packages #
#################
_java=$(command -v java)

if [[ "$_java" ]]; then
    java_ver=$("$_java" -version 2>&1 | awk -F '"' '/version/ {print $2}')
    if ! [[ "$java_ver" > "1.8" ]];
    then
        echo "Found java ${java_ver}, need 1.8 or greater."
        echo "Installation of java8 ..."
        brew cask install caskroom/versions/java8 >> ${LOGFILE}
        if [ $? -ne 0 ];
        then
            echo "Installation of java8 failed see ${LOGFILE} for details."
            exit 3
        fi
    fi
else
    echo "No java found." | tee -a ${LOGFILE}
    echo "Installation of java8 ..." | tee -a ${LOGFILE}
    brew cask install caskroom/versions/java8 >> ${LOGFILE} 2>&1
    if [ $? -ne 0 ];
    then
        echo "Installation of java8 failed see ${LOGFILE} for details."
        exit 3
    fi
fi

for item in "${brew_packages[@]}";
do
    pack=${item%%:*}
    cmd=${item#*:}
    brew_install "${pack}" "${cmd}" || exit 3
done


#################
# Perl packages #
#################

echo "Installation of Perl and packages"

if ! command_exists clustalw;
then
   echo "Create link to clustalw needed by Bioperl:clustalw" | tee -a ${LOGFILE}
   clustalw2="$(command -v clustalw2)"
   sudo ln -s ${clustalw2} ${PREFIX}bin/clustalw || exit 2
fi

if [[ -d ${crispr_perllib} ]];
then
    sudo mkdir -p ${crispr_perllib}
fi

for pack in ${perl_packages} ;
do
    if ! perl_pack_exists ${pack};
    then
        echo "Installation of ${pack}..." | tee -a ${LOGFILE}
        sudo cpanm -l ${crispr_perllib} ${pack} >> ${LOGFILE} 2>&1
        if [ $? -ne 0 ];
        then
            echo "Installation of ${pack} failed see ${LOGFILE} for details."
            exit 4
        fi
    fi
done

for pack in ${perl_packages_forced} ;
do
    if ! perl_pack_exists ${pack};
    then
        echo "Installation of ${pack}." | tee -a ${LOGFILE}
        sudo cpanm -l ${crispr_perllib} --force ${pack} >> ${LOGFILE} 2>&1
        if [ $? -ne 0 ];
        then
            echo "Installation of ${pack} failed see ${LOGFILE} for details."
            exit 4
        fi
    fi
done


##########
# vmatch #
##########
if ! command_exists vmatch2 ;
then
    echo "Installation of Vmatch..." | tee -a ${LOGFILE}
    if [ ! -f ${vmatch}.tar.gz ];
    then
        download vmatch ${vmatch_url} >> ${LOGFILE} 2>&1
    fi
    tar -zxf ${vmatch}.tar.gz >> ${LOGFILE} 2>&1
    gcc -Wall -Werror -fPIC -O3 -shared ${vmatch}/SELECT/sel392.c -o sel392v2.so >> ${LOGFILE} 2>&1

    test -d ${PREFIX}lib/ || mkdir -p ${PREFIX}lib/
    sudo install -m 0775 sel392v2.so ${PREFIX}lib/ >> ${LOGFILE} 2>&1 && \
    sudo install -m 0775 ${vmatch}/vmatch ${PREFIX}bin/vmatch2 >> ${LOGFILE} 2>&1 && \
    sudo install -m 0775 ${vmatch}/vsubseqselect ${PREFIX}bin/vsubseqselect2 >> ${LOGFILE} 2>&1 && \
    sudo install -m 0775 ${vmatch}/mkvtree ${PREFIX}bin/mkvtree2 >> ${LOGFILE} 2>&1 &&\
    if [ $? -ne 0 ];
        then
            echo "Installation of vmatch failed see ${LOGFILE} for details."
            exit 8
        fi
fi

###############
# MacSyFinder #
###############
if ! command_exists macsyfinder ;
then
    echo "Installation of MacSyFinder" | tee -a ${LOGFILE}
    if [ ! -f ${macsyfinder}.tar.gz ];
    then
        download macsyfinder ${macsyfinder_url} | tee -a ${LOGFILE}
    fi
    tar -xzf ${macsyfinder}.tar.gz
    cd ${macsyfinder}
    python2.7 setup.py build >> ${LOGFILE} 2>&1
    python2.7 setup.py test -vv >> ${LOGFILE} 2>&1
    if [ $? -ne 0 ];
        then
            echo "MacSyFinder tests failed see ${LOGFILE} for details."
            exit 3
        fi
    sudo python2.7 setup.py install >> ${LOGFILE} 2>&1
    if [ $? -ne 0 ];
    then
        echo "MacSyFinder installation failed see ${LOGFILE} for details."
        exit 5
    fi
    cd ${CURDIR}
fi

#############
# CasFinder #
#############

cas_data="${PREFIX}share/macsyfinder/"
cd ${CURDIR}

if [ ! -d ${cas_data} ];
then
    echo "Create directory for ${casfinder}" | tee -a ${LOGFILE}
    mkdir -p ${cas_data} >> ${LOGFILE} 2>&1
fi

if [ ! -d "${cas_data}/${casfinder}" ];
then
    echo "Copy ${casfinder} to ${cas_data}" | tee -a ${LOGFILE}
    sudo cp -pr "${casfinder}" "${cas_data}"
    if [ $? -ne 0 ];
    then
        echo "${casfinder} installation failed, see ${LOGFILE} for details."
        exit 5
    fi
else
    echo "${casfinder} found at ${cas_data}/${casfinder}." | tee -a ${LOGFILE}
    echo "Clean ${casfinder}." | tee -a ${LOGFILE}
    rm -Rf "${cas_data}/${casfinder}"
    # install cas profiles and definition packaged with CRISPRCasFinder
    echo "Install ${casfinder} from CRISPRCasFinder package." | tee -a ${LOGFILE}
    sudo cp -pr "${casfinder}" "${cas_data}"
fi

cd ${CURDIR}


#######################
# prokka dependencies #
#######################

###########
# tbl2asn #
###########

if ! command_exists macsyfinder ;
then
    echo "Installation of tbl2asn" | tee -a ${LOGFILE}
    tbl2asn_url="ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/mac.tbl2asn.gz"
    # curl does not work ;-(
    wget "${tbl2asn_url}" || echo "Cannot download ${tbl2asn_url} download it manually and run again this script" && exit 1
    gunzip mac.tbl2asn.gz
    sudo install -m 0755 mac.tbl2asn ${PREFIX}bin/tbl2asn
fi

##########
# prokka #
##########

if ! command_exists prokka ;
then
    echo "Installation of prokka" | tee -a ${LOGFILE}
     if [ ! -f ${prokka}.tar.gz ];
    then
        download macsyfinder "${prokka_url}" | tee -a ${LOGFILE}
    fi
    tar -xzf ${prokka}.tar.gz
    cd ${prokka}

    prokka_data="${PREFIX}share/prokka"
    prokka_db="${prokka_data}/db"
    test -d "${prokka_db}" || sudo mkdir -p "${prokka_db}"
    # copy database
    sudo cp -pr db/* ${prokka_db}

    # tell prokka where to find its tools and db once installed
    sed -i -e "s|my \$BINDIR.*|my \$BINDIR=\"${PREFIX}libexec/prokka\";|" bin/prokka || exit 6
    sed -i -e "s|my \$DBDIR.*|my \$DBDIR=\"${prokka_db}\";|" bin/prokka || exit 6


    for _bin in bin/*;
    do
        sudo install -m 0755 ${_bin} ${PREFIX}bin/
    done

    # install prokka binaries
    test -d ${PREFIX}/libexec/prokka || sudo mkdir -p ${PREFIX}libexec/prokka

    for p in binaries/darwin/*;
    do
        sudo install -m 0755 ${p} ${PREFIX}libexec/prokka || exit 6
    done
    # parallel is installed via packet manager
    sudo install -m 0755 binaries/common/minced ${PREFIX}libexec/prokka/ && \
    sudo install -m 0644 binaries/common/minced.jar ${PREFIX}libexec/prokka/ || exit 6

    # setup prokka db
    prokka_cmd="${PREFIX}bin/prokka"

    sudo ${prokka_cmd} --setupdb
    ${prokka_cmd} --version
    cd "${CURDIR}"
fi

###################
# CRISPRCasFinder #
###################

crispr_data=${PREFIX}/share/CRISPRCasFinder

patch CRISPRCasFinder.pl CRISPRCasFinder.patch && \
sudo install -m 0755 CRISPRCasFinder.pl ${PREFIX}/bin/CRISPRCasFinder && \
sudo install -m 0644 supplementary_files/crispr.css ${crispr_data} && \
sudo install -m 0644 supplementary_files/Repeat_List.csv ${crispr_data} && \
sudo install -m 0644 supplementary_files/CRISPR_crisprdb.csv ${crispr_data} && \
sudo install -m 0644 supplementary_files/repeatDirection.tsv ${crispr_data}

if [ $? -ne 0 ];
then
    echo "Cannot install CRISPRCasFinder see ${LOGFILE} for details."
    exit 9
fi

if test $?;
then
    echo "Installation finished."
fi

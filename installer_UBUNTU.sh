#!/bin/bash

# shell script allowing to install all CRISPRCasFinder.pl-v4.2's dependencies
#
# same version than CRISPRCasFinder.pl, here 4.2.17
# authors: David Couvin, Fabrice Leclerc, Claire Toffano-Nioche

#------------------------
function launchInstall {
# $1 $packageManagmentInstall
# $2 package name
# $3 $LOGFILE
    echo "Installation of $2" >> $3
    $1 $2 >> $3
}
#------------------------

#save the current directory:
CURDIR=`pwd` 

#create log file:
LOGFILE=$CURDIR/install.log
if [ -e $LOGFILE ]; then rm $LOGFILE ; fi

sudo chmod 755 .
sudo chmod 755 *

#create src & bin folders:
echo "create $CURDIR/src and $CURDIR/bin folders" >> $LOGFILE
if [ ! -d $CURDIR/src ];then mkdir $CURDIR/src; fi
if [ ! -d $CURDIR/bin ];then mkdir $CURDIR/bin; fi

#geting OS info: 
ostype=`echo $OSTYPE`
echo "$ostype" >> $LOGFILE

if [ ! "$ostype" = "linux-gnu" ]; then
    echo 'Sorry, install process only for OSTYPE linux-gnu (based on apt-get)'
else
    echo "apt-get update/upgrade" >> $LOGFILE
    sudo apt-get -y update >> $LOGFILE
    sudo apt-get -y upgrade >> $LOGFILE
    packageManagmentInstall='sudo apt-get -y install '
    distribution='Linux_x86_64'

    #important packages
    launchInstall "$packageManagmentInstall" "wget" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "curl" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "git" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "default-jre" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "python" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "parallel" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "cpanminus" "$LOGFILE"

    #"bioinfo" packages
    launchInstall "$packageManagmentInstall" "hmmer" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "emboss" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "emboss-lib" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "ncbi-blast+" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "bioperl" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "bioperl-run" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "libdatetime-perl" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "libxml-simple-perl" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "libdigest-md5-perl" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "clustalw" "$LOGFILE"

    echo "copy /usr/bin/clustalw to $CURDIR/bin/clustalw2" >> $LOGFILE
    sudo cp /usr/bin/clustalw $CURDIR/bin/clustalw2
    sudo cpanm Try::Tiny >> $LOGFILE
    sudo cpanm Test::Most >> $LOGFILE
    sudo cpanm JSON::Parse >> $LOGFILE
    sudo cpanm Class::Struct >> $LOGFILE
    sudo cpanm Bio::DB::Fasta >> $LOGFILE
    sudo cpanm File::Copy  >> $LOGFILE
    sudo cpanm Bio::Seq Bio::SeqIO >> $LOGFILE
    sudo cpanm Bio::Tools::Run::Alignment::Clustalw >> $LOGFILE 
    sudo cpanm Bio::Tools::Run::Alignment::Muscle >> $LOGFILE

    #install vmatch
    echo "Installation of Vmatch" >> $LOGFILE
    echo "change directory to $CURDIR/src" >> $LOGFILE
    cd $CURDIR/src
    wget http://vmatch.de/distributions/vmatch-2.3.0-${distribution}-64bit.tar.gz >> $LOGFILE
    tar -zxf vmatch-2.3.0-${distribution}-64bit.tar.gz >> $LOGFILE
    gcc -Wall -Werror -fPIC -O3 -shared vmatch-2.3.0-${distribution}-64bit/SELECT/sel392.c -o $CURDIR/sel392v2.so >> $LOGFILE
    echo "copy $CURDIR/src/vmatch-2.3.0-${distribution}-64bit/vmatch, mkvtree and vsubseqselect to $CURDIR/bin/" >> $LOGFILE
    sudo cp $CURDIR/src/vmatch-2.3.0-${distribution}-64bit/vmatch $CURDIR/bin/vmatch2
    sudo cp $CURDIR/src/vmatch-2.3.0-${distribution}-64bit/mkvtree $CURDIR/bin/mkvtree2
    sudo cp $CURDIR/src/vmatch-2.3.0-${distribution}-64bit/vsubseqselect $CURDIR/bin/vsubseqselect2
    echo "change directory to $CURDIR" >> $LOGFILE
    cd $CURDIR

    #install muscle
    launchInstall "$packageManagmentInstall" "muscle" "$LOGFILE"
 
    #install prodigal
    launchInstall "$packageManagmentInstall" "prodigal" "$LOGFILE"

    #install macsyfinder
    echo "Installation of MacSyFinder" >> $LOGFILE
    cd ${CURDIR}
    wget https://dl.bintray.com/gem-pasteur/MacSyFinder/macsyfinder-1.0.5.tar.gz >> $LOGFILE
    tar -xzf macsyfinder-1.0.5.tar.gz
    test -d bin ||  mkdir bin
    cd bin
    ln -s ../macsyfinder-1.0.5/bin/macsyfinder
    cd ${CURDIR}
    echo "add definition of MACSY_HOME (${CURDIR}/macsyfinder-1.0.5/) in .profile" >> $LOGFILE
    echo "export MACSY_HOME=${CURDIR}/macsyfinder-1.0.5/" >> $HOME/.profile

    echo "add bin folder ($CURDIR/bin) to the definition of PATH in $HOME/.profile" >> $LOGFILE
    echo "export PATH=${CURDIR}/bin:${PATH}" >> $HOME/.profile
    
    #set environment variables
    source $HOME/.profile

# if $OSTYPE
fi 


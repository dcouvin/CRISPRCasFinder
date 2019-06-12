#!/bin/bash

# shell script allowing to install all CRISPRCasFinder.pl-v4.2's dependencies on Fedora
# (could be potentially adapted for Red Hat)
#
# same version than CRISPRCasFinder.pl, here 4.2.19
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
    echo 'Sorry, install process only for OSTYPE linux-gnu (based on yum)'
else
    echo "yum update/upgrade" >> $LOGFILE
    #sudo yum -y update >> $LOGFILE
    #sudo yum -y upgrade >> $LOGFILE
    packageManagmentInstall='sudo yum -y install '
    distribution=''

    if $(uname -m | grep '64'); then
       distribution='Linux_x86_64-64bit'
       echo "$ostype distribution is: 64-bit" >> $LOGFILE
    else
       distribution='Linux_i386-32bit'
       echo "$ostype distribution is: 32-bit" >> $LOGFILE
    fi

    #important packages
    if [ ! -x "$(command -v wget)" ]; then
        launchInstall "$packageManagmentInstall" "wget" "$LOGFILE"
    fi
    if [ ! -x "$(command -v curl)" ]; then
        launchInstall "$packageManagmentInstall" "curl" "$LOGFILE"
    fi
    if [ ! -x "$(command -v git)" ]; then
        launchInstall "$packageManagmentInstall" "git" "$LOGFILE"
    fi
    if [ ! -x "$(command -v java)" ]; then
        launchInstall "$packageManagmentInstall" "java" "$LOGFILE"
    fi
    if [ ! -x "$(command -v python2.7)" ]; then
        launchInstall "$packageManagmentInstall" "python-unversioned-command" "$LOGFILE"
    fi
    if [ ! -x "$(command -v parallel)" ]; then
        launchInstall "$packageManagmentInstall" "parallel" "$LOGFILE"
    fi
    if [ ! -x "$(command -v cpanm)" ]; then
        launchInstall "$packageManagmentInstall" "perl-App-cpanminus" "$LOGFILE"
    fi
    if [ ! -x "$(command -v cpan)" ]; then
        launchInstall "$packageManagmentInstall" "perl-CPAN" "$LOGFILE"
    fi
    launchInstall "$packageManagmentInstall" "perl-Module-Build" "$LOGFILE"
  
    #"bioinfo" packages
    if [ ! -x "$(command -v hmmsearch)" ]; then
        launchInstall "$packageManagmentInstall" "hmmer" "$LOGFILE"
        #launchInstall "$packageManagmentInstall" "hmmer-devel" "$LOGFILE"   # may not be necessary
    fi  
    if [ ! -x "$(command -v fuzznuc)" ] && [ ! -x "$(command -v needle)" ]
    then
        launchInstall "$packageManagmentInstall" "EMBOSS" "$LOGFILE"  
        #launchInstall "$packageManagmentInstall" "EMBOSS-devel" "$LOGFILE"  # may not be necessary
        launchInstall "$packageManagmentInstall" "EMBOSS-libs" "$LOGFILE"
    fi
	
    launchInstall "$packageManagmentInstall" "ncbi-blast+" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "perl-Time-Piece" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "perl-XML-Simple" "$LOGFILE"
    launchInstall "$packageManagmentInstall" "perl-Digest-MD5" "$LOGFILE"
    
	if [ ! -x "$(command -v clustalw)" ]; then
        launchInstall "$packageManagmentInstall" "clustalw" "$LOGFILE"  # should work for CentOS
    fi
    if [ ! -x "$(command -v gcc)" ]; then
        launchInstall "$packageManagmentInstall" "gcc-c++" "$LOGFILE"
    fi
    if [ ! -x "$(command -v make)" ]; then
        launchInstall "$packageManagmentInstall" "make" "$LOGFILE"
    fi

    #clustalw
    #test if '/usr/bin/clustalw' file exists, otherwise install clustalw manually
    CLUSTALWFILE=/usr/bin/clustalw
    if [ -f "$CLUSTALWFILE" ]
    then
        echo "copy /usr/bin/clustalw to $CURDIR/bin/clustalw2" >> $LOGFILE
        sudo cp /usr/bin/clustalw $CURDIR/bin/clustalw2
    else
        #manual installation of clustalw
        wget http://www.clustal.org/download/current/clustalw-2.1.tar.gz >> $LOGFILE
        tar -xzvf clustalw-2.1.tar.gz >> $LOGFILE
        cd clustalw-2.1
		./configure >> $LOGFILE
		sudo make >> $LOGFILE
        sudo make install >> $LOGFILE
        sudo cp src/clustalw2 $CURDIR/bin/clustalw2
		sudo cp src/clustalw2 /usr/bin/clustalw
        cd $CURDIR
    fi

    #cpanm
    sudo cpanm Bio::Perl >> $LOGFILE
    sudo cpanm Bio::FeatureIO >> $LOGFILE
    sudo cpanm Try::Tiny >> $LOGFILE
    sudo cpanm Test::Most >> $LOGFILE
    sudo cpanm JSON::Parse >> $LOGFILE
    sudo cpanm Class::Struct >> $LOGFILE
    sudo cpanm Bio::DB::Fasta >> $LOGFILE
    sudo cpanm File::Copy  >> $LOGFILE
    sudo cpanm Bio::Seq Bio::SeqIO >> $LOGFILE
    sudo cpanm --force Bio::Tools::Run::Alignment::Clustalw >> $LOGFILE 
    sudo cpanm --force Bio::Tools::Run::Alignment::Muscle >> $LOGFILE
    sudo cpanm Date::Calc >> $LOGFILE

    #BioPerl-Run
    # wget https://cpan.metacpan.org/authors/id/C/CJ/CJFIELDS/BioPerl-Run-1.007003.tar.gz 
    # most recent version of BioPerl-Run will not work. We should use release-1-7-0
    wget https://github.com/bioperl/bioperl-run/archive/release-1-7-0.tar.gz >> $LOGFILE
    tar -xzvf release-1-7-0.tar.gz >> $LOGFILE
    cd bioperl-run-release-1-7-0
    perl Build.PL
    sudo ./Build test
	sudo ./Build installdeps
    sudo ./Build install
    cd $CURDIR
	
    #install vmatch
    if [ ! -x "$(command -v vmatch2)" ] && [ ! -x "$(command -v mkvtree2)" ] && [ ! -x "$(command -v vsubseqselect2)" ]
    then
      echo "Installation of Vmatch" >> $LOGFILE
      echo "change directory to $CURDIR/src" >> $LOGFILE
      cd $CURDIR/src
      wget http://vmatch.de/distributions/vmatch-2.3.0-${distribution}.tar.gz >> $LOGFILE
      tar -zxf vmatch-2.3.0-${distribution}.tar.gz >> $LOGFILE
      gcc -Wall -Werror -fPIC -O3 -shared vmatch-2.3.0-${distribution}/SELECT/sel392.c -o $CURDIR/sel392v2.so >> $LOGFILE
      echo "copy $CURDIR/src/vmatch-2.3.0-${distribution}/vmatch, mkvtree and vsubseqselect to $CURDIR/bin/" >> $LOGFILE
      sudo cp $CURDIR/src/vmatch-2.3.0-${distribution}/vmatch $CURDIR/bin/vmatch2
      sudo cp $CURDIR/src/vmatch-2.3.0-${distribution}/mkvtree $CURDIR/bin/mkvtree2
      sudo cp $CURDIR/src/vmatch-2.3.0-${distribution}/vsubseqselect $CURDIR/bin/vsubseqselect2
      echo "change directory to $CURDIR" >> $LOGFILE
      cd $CURDIR
    else
        echo "Vmatch2, mkvtree2, and vsubseqselect2 are already installed." >> $LOGFILE
    fi

    #install muscle
    if [ ! -x "$(command -v muscle)" ]; then
        echo "Installation of Muscle" >> $LOGFILE
        launchInstall "$packageManagmentInstall" "muscle" "$LOGFILE" # may not work properly
    else
        echo "Muscle is already installed." >> $LOGFILE
    fi
    
    if [ ! -x "$(command -v muscle)" ]; then
      echo "Installation of Muscle (second trial)" >> $LOGFILE
      if [ $distribution == 'Linux_x86_64-64bit' ]
      then
        wget http://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz >> $LOGFILE
        tar -xzvf muscle3.8.31_i86linux64.tar.gz >> $LOGFILE
        sudo chmod +x muscle3.8.31_i86linux64
        sudo cp muscle3.8.31_i86linux64 /usr/bin/muscle
      else 
        wget http://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux32.tar.gz >> $LOGFILE
        tar -xzvf muscle3.8.31_i86linux32.tar.gz >> $LOGFILE
        sudo chmod +x muscle3.8.31_i86linux32
        sudo cp muscle3.8.31_i86linux32 /usr/bin/muscle
      fi
    else
        echo "Muscle is already installed." >> $LOGFILE
    fi

    #install prodigal
    if [ ! -x "$(command -v prodigal)" ] 
    then
        echo "Installation of Prodigal" >> $LOGFILE
        wget https://github.com/hyattpd/Prodigal/archive/GoogleImport.zip >> $LOGFILE
        unzip GoogleImport.zip
        cd Prodigal-GoogleImport
        sudo make install
    else
        echo "Prodigal is already installed." >> $LOGFILE
    fi

    #install macsyfinder
    if [ ! -x "$(command -v macsyfinder)" ] 
    then
      echo "Installation of MacSyFinder" >> $LOGFILE
      cd ${CURDIR}
      wget https://dl.bintray.com/gem-pasteur/MacSyFinder/macsyfinder-1.0.5.tar.gz >> $LOGFILE
      tar -xzf macsyfinder-1.0.5.tar.gz
      test -d bin ||  mkdir bin
      cd bin
      ln -s ../macsyfinder-1.0.5/bin/macsyfinder
      cd ${CURDIR}
      echo "add definition of MACSY_HOME (${CURDIR}/macsyfinder-1.0.5/) in .bashrc" >> $LOGFILE
      echo "export MACSY_HOME=${CURDIR}/macsyfinder-1.0.5/" >> $HOME/.bashrc
      echo "add bin folder ($CURDIR/bin) to the definition of PATH in $HOME/.bashrc" >> $LOGFILE
      echo "export PATH=${CURDIR}/bin:${PATH}" >> $HOME/.bashrc
    else
      echo "MacSyFinder is already installed." >> $LOGFILE
    fi

    #set environment variables
    #source $HOME/.bashrc #this command must be typed directly by user
    
    echo "Installation done."
	echo "You may need to reinstall some Perl's modules (with command: sudo cpanm ...)"
	echo "Please type command 'source ~/.bashrc' in order to reload your .bashrc file, then run 'perl CRISPRCasFinder.pl -v' to see if everything is correct."
# if $OSTYPE
fi 

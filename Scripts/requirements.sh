#!/bin/bash

if [[ "$OSTYPE" == linux-gnu* ]]
then
    #Download anaconda file for Linux
    wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
    #Checking the integrity of the file
    sha256sum Anaconda3-2020.07-Linux-x86_64.sh
    #Running the .sh script
    bash Anaconda3-2020.07-Linux-x86_64.sh
    #Compiling from source
    source ~/.bashrc
    rm Anaconda3-2020.07-Linux-x86_64.sh
    #Installing R
    sudo apt update
    sudo apt -y upgrade
    sudo apt -y install r-base
    #Installing R packages
    sudo apt-get install libcurl4-openssl-dev libxml2-dev
    sudo apt-get install libssl-dev
    sudo add-apt-repository -y ppa:cran/imagemagick
    sudo apt-get update
    sudo apt-get install -y libmagick++-dev
    sudo apt-get install gedit
    #Install Python packages
    python3 -m pip install pysamstats
    python3 -m pip install pysams
fi


if [[ "$OSTYPE" = darwin* ]]
then
    #Download Homebrew 
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
    #Download anaconda and miniconda 
    brew install --cask anaconda
    brew cask install miniconda
    #install wget
    brew install wget
    #Install R
    brew install r
    #Install Python packages
    pip install pysamstats
    pip install pysam
fi


#Using conda to install bowtie2, samtools, bedtools, picard.
conda install -c bioconda bowtie2=2.3.4.2
conda install -c bioconda samtools=1.9
conda install -c bioconda bedtools=2.28.0
conda install -c bioconda picard=2.22.1

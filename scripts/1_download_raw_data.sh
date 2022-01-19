#!/bin/bash
#SBATCH --job-name=downlaod_raw_data
#SBATCH --output=download_raw_data_log.txt
#SBATCH --time=3:30:00

#load required models
module load vital-it/latest
module load UHTS/Analysis/sratoolkit/2.10.7 
#Please run "vdb-config --interactive" command to configure sratoolkit before using it.

#change directory to raw data folder
cd ../data/raw

#download samples using SRA prefatch
prefetch --option-file samplesName.txt 

#convert sra files to fastq using fastq dump. Please use --gzip	or --bzip2 to compress the output fastq files
fastq-dump */*.sra --gzip

#delete sra file folders
rm -r */

#rename files

mv SRR9596295.fastq.gz somata_1.fastq.gz
mv SRR9596296.fastq.gz somata_2.fastq.gz
mv SRR9596300.fastq.gz somata_3.fastq.gz
mv SRR9596310.fastq.gz neuropil_1.fastq.gz
mv SRR9596303.fastq.gz neuropil_2.fastq.gz
mv SRR9596304.fastq.gz neuropil_3.fastq.gz





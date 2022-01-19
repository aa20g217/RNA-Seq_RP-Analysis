#!/bin/bash

#SBATCH --job-name=clipp_trimm_fastqc_multiqc
#SBATCH --output=clipp_trimm_fastqc_multiqc_log.txt
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:30:00
#SBATCH --mem-per-cpu=32G

#load required models
module load vital-it/latest
module load UHTS/Quality_control/cutadapt/2.5
module load UHTS/Quality_control/fastqc/0.11.9
module load UHTS/Analysis/MultiQC/1.8

#cd to raw data file folder
cd ../data/raw

###################Â
#iterate through each raw fastq file and remove adapters
################### 

#Clip "AGATCGGAAGAGCACACGTCTGAA" sequence from the 3' end.
#Cut 2 nt from the 5' end


for filename in *.gz;do

cutadapt ${filename} -j 32 -a AGATCGGAAGAGCACACGTCTGAA -o "${filename%.fastq.gz}_clipped.fastq.gz" \
--discard-untrimmed  --cut 2 -q 25 --minimum-length 22 --overlap 3 -e 0.2

done

################### 
#iterate through each clipped fastq file and trim them
################### 
#cut 10 nt from the 3' end


for filename in *_clipped.fastq.gz;do

cutadapt ${filename} -j 32 -o "${filename%.fastq.gz}_trimmed.fastq.gz" -q 25 \
--cut -10 --minimum-length 22 

done


################### 
#run fastqc and multiQC for clipped and trimmed files
################### 

mkdir ../QC						#make dir to store fastQC result for clipped files
fastqc *_clipped.fastq.gz -o ../QC/ -t 32 				# runfastqc


mkdir ../QC                                              #make dir to store fastQC result for trimmed files
fastqc *trimmed* -o ../QC/ -t 32                          # runfastqc


multiqc ../QC/*fastqc.zip \
	-o ../QC/ --title "raw_clipped_trimmed_reads_fastqc_multiqc_report"          # run multiqc to create a combine report

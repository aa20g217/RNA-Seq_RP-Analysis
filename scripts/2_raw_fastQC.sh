#!/bin/bash
#
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc_log.txt

#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1

#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=2

#SBATCH --time=02:00:00

#load required models
module load vital-it/latest
module load UHTS/Quality_control/fastqc/0.11.9 

mkdir ../data/QC 					#create a dir to save fastqc report

#run fastqc for every fastq file in raw data folder
fastqc -o ../data/QC/ -t 6 ../data/raw/*.fastq.gz

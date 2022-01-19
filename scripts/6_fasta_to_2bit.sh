#!/bin/bash

#SBATCH --job-name=fato2bit
#SBATCH --output=fato2bitc_log.txt
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:30:00
#SBATCH --mem-per-cpu=32G

#load required models
module load vital-it/latest
module load SequenceAnalysis/blat/36 

#cd to ref genomme and convert it to 2 bit
cd ../data/refGenome/

faToTwoBit Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa Rattus_norvegicus.Rnor_6.0.dna.toplevel.2bit


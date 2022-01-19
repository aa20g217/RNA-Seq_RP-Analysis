#!/bin/bash

#SBATCH --job-name=mapping
#SBATCH --output=mappingc_log.txt
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=05:30:00
#SBATCH --mem-per-cpu=32G

#load required models
module load vital-it/latest
module load UHTS/Aligner/bowtie/1.2.0
module add UHTS/Analysis/samtools/1.10

#cd to raw data file folder
cd ../data/raw


###################?
#iterate through each trimmed fastq file and keep only those read that does
# not mapped to undesired reads using -un feature from bowtie
###################

#unzipp fastq files
for filename in *_trm.fastq.gz;do

gunzip ${filename}
done


echo "***********undesired rna*****************"
for filename in *_trm.fastq;do

echo ${filename}

bowtie -S -t -p 12 ../index/Rnor_r_sno_sn_t_RNA  ${filename} \
--un "${filename%.fastq}_desired_reads.fastq" \
2> "${filename%.fastq}_desired_reads_log.txt" > /dev/null

done


###################?
#iterate through each desire read fastq file and 
# map them to ref genome
###################

#make directory to save alligned reads
mkdir ../allignedReads

echo "***********desired rna*****************"

for filename in *_desired_reads.fastq;do

echo ${filename}

bowtie -S -t -p 12 -v 1 -m 1 --best --strata \
../index/Rattus_norvegicus.Rnor_6.0.dna.toplevel -q ${filename} | \

samtools view -h -F 4 -b > ../allignedReads/"${filename%.fastq}_Rnor_6.bam";

done


###################?
#iterate through each bam file and sort them
###################

cd ../allignedReads

for filename in *.bam;do
samtools sort -@ 12 ${filename} -o "${filename%.bam}_sorted.bam" 

done


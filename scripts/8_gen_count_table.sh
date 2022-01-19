#!/bin/bash

#SBATCH --job-name=gen_count_table
#SBATCH --output=gen_count_table.txt
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:30:00
#SBATCH --mem-per-cpu=32G

#load required models
module load vital-it/latest
module add UHTS/Analysis/subread/2.0.1
module load UHTS/Analysis/MultiQC/1.8

#creat  countsTable dir
mkdir ../data/countsTable

#setup paths
cd ..
echo ${PWD}

# Count reads on CDS
featureCounts -T 4 -t CDS -g gene_id -a data/refGenome/Rattus_norvegicus.Rnor_6.0.104.gtf -o data/countsTable/CDS_counts_rawfile.txt data/allignedReads/*_sorted.bam

# process count file
cut -f 1,7-12 data/countsTable/CDS_counts_rawfile.txt > data/countsTable/CDS_counts_processed.txt



# Extract reads mapped to different biotypes
featureCounts -T 4 -t exon -g gene_biotype -a data/refGenome/Rattus_norvegicus.Rnor_6.0.104.gtf -o data/countsTable/biotype_counts_rawfile.txt data/allignedReads/*_sorted.bam

# Extract Biotype and Sample columns
cut -f 1,7-12 data/countsTable/biotype_counts_rawfile.txt > data/countsTable/biotype_counts_processed.txt

multiqc CDS_counts_rawfile.txt.summary --title "quantification_cds"
multiqc biotype_counts_rawfile.txt.summary --title "quantification_bio"


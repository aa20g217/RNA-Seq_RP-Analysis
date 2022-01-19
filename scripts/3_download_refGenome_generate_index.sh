#!/bin/bash
#
#SBATCH --job-name=download_refGenome_generate_index
#SBATCH --output=download_refGenome_generate_index_log.txt
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:30:00
#SBATCH --mem-per-cpu=32G

#load required models
module load vital-it/latest
module load UHTS/Aligner/bowtie/1.2.0 

#combine multiple undesired rna files into one
cat ../data/undesiredRNA/*.txt > ../data/undesiredRNA/Rnor_r_sno_sn_t_RNA.fa

#create dir to save ref genome
mkdir ../data/refGenome

#change directory to reference genome
cd ../data/refGenome

#download ref genome and annotaton file
wget http://ftp.ensembl.org/pub/release-104/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
wget http://ftp.ensembl.org/pub/release-104/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.104.gtf.gz

#unzip downloaded file and remove zipped file
gunzip *.gz 
rm *.gz

#create	dir to save genome index 
mkdir ../index

#create index for undesired RNA
bowtie-build  ../undesiredRNA/Rnor_r_sno_sn_t_RNA.fa ../index/Rnor_r_sno_sn_t_RNA -p 16

#create	index for genome
bowtie-build  ../refGenome/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa ../index/Rattus_norvegicus.Rnor_6.0.dna.toplevel -p 16

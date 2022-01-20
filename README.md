## Scripts Required to Analyse Ribosomal Profiling Data
### Requirements
Used Tools and version:
```
1) vital-it : latest
2) sratoolkit : 2.10.7
3) FastQC : 0.11.9
4) MultiQC : 1.8
5) Bowtie: 1.2.0 
6) Cutadapt: 2.5
7) Python: 3.5.2
8) Samtools: 1.10	
9) BLAT: 36
10) R:	4.0.5
11) RiboseQC: 0.99.0 (please see sessionInfo file for dependencies)
12) Subread: 2.0.1
``` 

### Reference Genome
1) Genome: http://ftp.ensembl.org/pub/release-104/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
2) Annotation: http://ftp.ensembl.org/pub/release-104/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.104.gtf.gz

### How to Run 
All the required scripts are available in ```scripts folder```. Please run them in the given order. All the results will be saved in ```data folder``` and logs will be available in scripts folder.

### Notes:
1) Please run ```vdb-config --interactive``` command to configure sratoolkit before using it.
2) Please refer to session info files in the script folder for all the dependencies of RiboSeQC and DESeq2.

# Load the package
library("RiboseQC")

setwd("/Users/akshay/Desktop/PhD/Course-work/RNA-seq/lectures/ribosome profiling module/refGenome")

#prepare annotation file
prepare_annotation_files(annotation_directory = "../genome_for_RiboseQC/",
                         gtf_file = "Rattus_norvegicus.Rnor_6.0.104.gtf",
                         twobit_file = "Rattus_norvegicus.Rnor_6.0.dna.toplevel.2bit",
                         scientific_name = "Rattus.norvegicus",
                         annotation_name = "Rnor_6",
                         export_bed_tables_TxDb = F)



#create a list of sorted BAM fiiles
setwd("/Users/akshay/Desktop/PhD/Course-work/RNA-seq/lectures/ribosome profiling module/sortedBam")
genome_bam =list.files()


# create QC report 
RiboseQC_analysis(annotation_file ="../genome_for_RiboseQC/Rattus_norvegicus.Rnor_6.0.104.gtf_Rannot",
                  bam_files = genome_bam,
                  report_file = "Rnor_Riboseq_QC.html",
                  sample_names = c("Neuropil_Poly_1", "Neuropil_Poly_2", "Neuropil_Poly_3",
                                   "Somata_Poly_1", "Somata_Poly_2", "Somata_Poly_3"),
                  dest_names = c("Neuropil_Poly_1", "Neuropil_Poly_2", "Neuropil_Poly_3",
                                 "Somata_Poly_1", "Somata_Poly_2", "Somata_Poly_3"),
                  write_tmp_files = F)


setwd("/Users/akshay/Desktop/PhD/Course-work/RNA-seq/lectures/ribosome profiling module/RiboseQC/")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

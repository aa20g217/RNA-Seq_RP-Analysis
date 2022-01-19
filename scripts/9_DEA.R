library(DESeq2)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(BiocParallel)
library(org.Rn.eg.db) 
register(MulticoreParam(4)) 

# Setting up color profiles from colorbrewer
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mycols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

sample_1 <- "neuropil"
sample_2 <- "somata"
no_of_reps <- 3

sample_column <- c(rep(sample_1, no_of_reps),
                   rep(sample_2, no_of_reps))

run_column <- c(paste(sample_1, "1", sep = "_"),
                paste(sample_1, "2", sep = "_"),
                paste(sample_1, "3", sep = "_"),
                paste(sample_2, "1", sep = "_"),
                paste(sample_2, "2", sep = "_"),
                paste(sample_2, "3", sep = "_"))

rep_column <- c("A", "B", "C",
                "A", "B", "C")

samples_df <- data.frame(sample_column,
                         run_column,
                         rep_column)

colnames(samples_df) <- c("sample", "run", "rep")

samples_df$condition <- factor(rep(c(sample_1, sample_2), each = no_of_reps))

rownames(samples_df) <- samples_df$run

# Load count data
setwd("/Users/akshay/Desktop/PhD/Course-work/RNA-seq/lectures/ribosome profiling module/countsTable/")
featurecount_data <- read.table("CDS_counts_processed.txt", header = TRUE, row.names = 1)

# Change colnames
# Make sure the column order in featurecount_data matches samples_df !!!
colnames(featurecount_data) <- rownames(samples_df)

#Import as DESeqDataSet (dds)
dds <- DESeqDataSetFromMatrix(countData = featurecount_data,
                              colData = samples_df,
                              design = ~ condition)

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Factor levels
# WT columns should be first followed by KO/treatment
dds$condition <- factor(c(rep(sample_1, no_of_reps),
                          rep(sample_2, no_of_reps)),
                        levels = c(sample_1,
                                   sample_2))

# Differential expression analysis
dds <- DESeq(dds)
colData(dds)


setwd("/Users/akshay/Desktop/PhD/Course-work/RNA-seq/lectures/ribosome profiling module/DEA/")

################################################################################
################################################################################
# QC
################################################################################
################################################################################

# Log transformation for data quality assessment
rld <- rlog(dds, blind = FALSE)

# Sample distance matrix
sampleDists <- as.matrix(dist(t(assay(rld))))
pdf("QC_sample_distance_matrix_CDS.pdf")
heatmap.2(as.matrix(sampleDists),
          key = T,
          trace = "none",
          col = colorpanel(100, "#2b8cbe", "#7bccc4", "white"),
          ColSideColors = mycols[dds$condition],
          RowSideColors = mycols[dds$condition],
          margin = c(10, 10), main = "Sample Distance Matrix")
dev.off()

# Count matrix heatmap
select <- order(rowMeans(counts(dds,normalized = TRUE)))
df <- as.data.frame(colData(dds)[ , c("condition","rep")])

pdf("QC_count_matrix_CDS.pdf")
pheatmap(assay(rld)[select,],
         cluster_rows = FALSE,
         show_rownames=FALSE,
         cluster_cols = FALSE,
         annotation_col = df)
dev.off()

# PCA plot

pdf("QC_PCA_CDS.pdf")
pcaData <- plotPCA(rld, intgroup = c("condition", "rep"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = rep)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()
################################################################################
################################################################################

################################################################################
################################################################################

# Using contrast so that the following code can be scaled with increasing number of samples
res <- results(dds, contrast = c("condition", sample_1, sample_2), alpha = 0.05)
as.data.frame(mcols(res, use.names = T))[,"description"]

# Adding gene names using org.Rn.eg.db
# Source: http://bioconductor.org/help/course-materials/2015/LearnBioconductorFeb2015/B02.1.1_RNASeqLab.html
# Also: https://support.bioconductor.org/p/66288/
# This function takes a list of IDs as first argument and their key type as the second argument.
# The third argument is the key type we want to convert to, the fourth is the AnnotationDb object to use.
# Finally, the last argument specifies what to do if one source ID maps to several target IDs:
# should the function return an NA or simply the first of the multiple IDs

convertIDs <- function( ids, from, to, db, ifMultiple = c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

# # Check columns in the database that you want to add:
columns(org.Rn.eg.db)

# Actual adding of the column
res$GeneID <- row.names(res)
res$gene_symbol <- convertIDs(row.names(res), "ENSEMBL", "SYMBOL", org.Rn.eg.db)

summary(res)

res_df <- as.data.frame(res)

# Which genes are translationally deregulated:
res_df$regulation_level <- ifelse((res_df$log2FoldChange > 0.5 & res_df$padj < 0.05), "Upregulated",
                                  ifelse((res_df$log2FoldChange < - 0.5 & res_df$padj < 0.05),
                                         "Downregulated", "Unchanged"))

write.table(res_df,
            file = "DESeq2_res.csv",
            sep = ",",
            row.names = F,
            col.names = T,
            quote = F)




# Plot

pdf("Volcano_plot.pdf", width = 6, height = 6)
ggplot(res_df,
       aes(x = -log10(padj),
           y = log2FoldChange,
           color = regulation_level)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("#404788FF", "#73D055FF", "#999999")) +
  xlab("-log10(adjusted p-value)") +
  ylab("Log2 fold change") +
  labs(color = "Regulation level") +
  theme_bw()+coord_flip()+ggtitle("Neuropil vs Somata", subtitle ="Threshold -> FC: ±0.5, adjPval: 0.05")
dev.off()



################################################################################
################################################################################
# GO analysis
################################################################################
################################################################################
library(clusterProfiler)
library(topGO)

sample_name = "Neuropil_Poly_vs_Somata_Poly_only_padj"

df <- res_df

rownames(df) <- df$GeneID
df <- df[order(df$padj), ]

# Define upregulated and downregulated genes based on padj value

genes_up <- which(df$padj < 0.05 & df$log2FoldChange > 0)
genes_down <- which(df$padj < 0.05 & df$log2FoldChange < 0)

all_genes_names <- rownames(df)

genes_up <- rownames(df)[genes_up]
genes_down <- rownames(df)[genes_down]

genelist_up <- factor(as.integer(all_genes_names %in% genes_up))
names(genelist_up) <- all_genes_names

genelist_down <- factor(as.integer(all_genes_names %in% genes_down))
names(genelist_down) <- all_genes_names

allGO2genes <- annFUN.org(whichOnto = "ALL",
                          feasibleGenes = NULL,
                          mapping = "org.Rn.eg.db",
                          ID = "ensembl")

GOdata_up_bp <- new("topGOdata",
                    ontology = "BP",
                    allGenes = genelist_up,
                    annot = annFUN.GO2genes,
                    GO2genes = allGO2genes,
                    nodeSize = 10)

GOdata_up_mf <- new("topGOdata",
                    ontology = "MF",
                    allGenes = genelist_up,
                    annot = annFUN.GO2genes,
                    GO2genes = allGO2genes,
                    nodeSize = 10)

GOdata_up_cc <- new("topGOdata",
                    ontology = "CC",
                    allGenes = genelist_up,
                    annot = annFUN.GO2genes, 
                    GO2genes = allGO2genes,
                    nodeSize = 10)

GOdata_down_bp <- new("topGOdata",
                      ontology = "BP",
                      allGenes = genelist_down,
                      annot = annFUN.GO2genes,
                      GO2genes = allGO2genes,
                      nodeSize = 10)

GOdata_down_mf <- new("topGOdata",
                      ontology = "MF",
                      allGenes = genelist_down,
                      annot = annFUN.GO2genes,
                      GO2genes = allGO2genes,
                      nodeSize = 10)

GOdata_down_cc <- new("topGOdata",
                      ontology = "CC",
                      allGenes = genelist_down,
                      annot = annFUN.GO2genes,
                      GO2genes = allGO2genes,
                      nodeSize = 10)


resultFis_up_bp <- runTest(GOdata_up_bp, statistic = "fisher")
resultFis_up_mf <- runTest(GOdata_up_mf, statistic = "fisher")
resultFis_up_cc <- runTest(GOdata_up_cc, statistic = "fisher")
resultFis_down_bp <- runTest(GOdata_down_bp, statistic = "fisher")
resultFis_down_mf <- runTest(GOdata_down_mf, statistic = "fisher")
resultFis_down_cc <- runTest(GOdata_down_cc, statistic = "fisher")

parse_tables <- function(GO_data, statistics)
{
  goEnrichment <- GenTable(GO_data, weightFisher = statistics, topNodes = 20)
  sub("< ", "", goEnrichment$weightFisher)
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep = ", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels = rev(goEnrichment$Term))
  goEnrichment$weightFisher <- as.numeric(sub("< ", "", goEnrichment$weightFisher))  
  goEnrichment
}

GOres_up_bp <- parse_tables(GOdata_up_bp, resultFis_up_bp)
GOres_up_mf <- parse_tables(GOdata_up_mf, resultFis_up_mf)
GOres_up_cc <- parse_tables(GOdata_up_cc, resultFis_up_cc)

GOres_down_bp <- parse_tables(GOdata_down_bp, resultFis_down_bp)
GOres_down_mf <- parse_tables(GOdata_down_mf, resultFis_down_mf)
GOres_down_cc <- parse_tables(GOdata_down_cc, resultFis_down_cc)


plot_GO <- function(GO_data, Ontology, Regulation, use_color) {
  GO_data$log_weightFisher <- (- log10(as.numeric(GO_data$weightFisher)))
  ggplot(GO_data, 
         aes(x = GO_data$log_weightFisher,
             y = GO_data$Term)) +
    geom_segment(aes(x = 0,
                     xend = GO_data$log_weightFisher,
                     y = GO_data$Term,
                     yend = GO_data$Term),
                 colour = use_color)  +
    geom_point(aes(size = GO_data$Significant),
               colour = use_color) +
    scale_size_area(name = "Gene counts") +
    xlab("Enrichment (- log10 Pvalue)") +
    ylab(Ontology) +
    ggtitle(Regulation) +
    scale_x_continuous() +
    theme_bw() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"))
}

plot_up_BP <- plot_GO(GOres_up_bp, "Biological Proccess", "TopGO Up (fisher's exact test)", "#404788FF")
plot_up_MF <- plot_GO(GOres_up_mf, "Molecular Function", "TopGO Up (fisher's exact test)", "#404788FF")
plot_up_CC <- plot_GO(GOres_up_cc, "Cellular Component", "TopGO Up (fisher's exact test)", "#404788FF")

# grid.arrange(grobs = list(p1,p2,p3))

plot_down_BP <- plot_GO(GOres_down_bp, "Biological Proccess", "TopGO Down (fisher's exact test)", "#73D055FF")
plot_down_MF <- plot_GO(GOres_down_mf, "Molecular Function", "TopGO Down (fisher's exact test)", "#73D055FF")
plot_down_CC <- plot_GO(GOres_down_cc, "Cellular Component", "TopGO Down (fisher's exact test)", "#73D055FF")

# grid.arrange(grobs = list(p1,p2,p3))

pdf(paste(sample_name, "Biological_Proccess_TopGO_Up_fisher.pdf", sep = "_"))
plot_up_BP
dev.off()

pdf(paste(sample_name, "Molecular_Function_TopGO_Up_fisher.pdf", sep = "_"))
plot_up_MF
dev.off()

pdf(paste(sample_name, "Cellular_Component_TopGO_Up_fisher.pdf", sep = "_"))
plot_up_CC
dev.off()

pdf(paste(sample_name, "Biological_Proccess_TopGO_down_fisher.pdf", sep = "_"))
plot_down_BP
dev.off()

pdf(paste(sample_name, "Molecular_Function_TopGO_down_fisher.pdf", sep = "_"))
plot_down_MF
dev.off()

pdf(paste(sample_name, "Cellular_Component_TopGO_down_fisher.pdf", sep = "_"))
plot_down_CC
dev.off()

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

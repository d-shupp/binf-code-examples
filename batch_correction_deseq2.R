# Load libraries
library(DESeq2)
library(sva)
library(PCAtools)
library(ggplot2)
library(patchwork)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(pheatmap)

# Set work directory
setwd("~/Documents/GMU/Research_Jafri/tcell_model/work")

# Load counts and metadata
txi <- readRDS("txi.rds")
counts_raw <- txi$counts
colData <- read.delim("samples.tsv")

# DESeq2 function
run_DESeq2 <- function(type, counts, colData) {
  # Create DESeq2 data set (dds)
  if (type == "txi") {
    ddsTxi <- DESeqDataSetFromTximport(counts, colData = colData, design = ~state)
    # Filter out low expression genes (less than 10 counts)
    keep <- rowSums(counts(ddsTxi)) >= 10
    ddsTxi <- ddsTxi[keep, ]
    # Set factor levels
    ddsTxi$state <- relevel(ddsTxi$state, ref = "Naive")
    dds <- DESeq(ddsTxi)
  } else if (type == "mat") {
    ddsMat <- DESeqDataSetFromMatrix(counts, colData = colData, design = ~state)
    # Filter out low expression genes (less than 10 counts)
    keep <- rowSums(counts(ddsMat)) >= 10
    ddsMat <- ddsMat[keep, ]
    # Set factor levels
    ddsMat$state <- relevel(ddsMat$state, ref = "Naive")
    dds <- DESeq(ddsMat)
  }
  return(dds)
}


##### Main #####


# Batch correction
counts_corrected <- ComBat_seq(counts_raw,
                               batch = colData$batch,
                               group = colData$state)
counts_corrected_rounded <- round(counts_corrected)

# Run analysis on raw counts
dds_raw <- run_DESeq2("txi", txi, colData)
# Normalization
vsd_raw <- assay(vst(dds_raw, blind = FALSE))

# Run analysis on batch-corrected counts
dds_corrected <- run_DESeq2("mat", counts_corrected_rounded, colData)
# Normalization
vsd_corrected <- assay(vst(dds_corrected, blind = FALSE))
vsd_corrected_ml <- assay(vst(dds_corrected, blind = TRUE))
# Save normalized data for ML pipeline
res <- results(dds_corrected)
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
top_genes <- head(sig_genes[order(sig_genes$padj), ], 1000)
top_genes_list <- rownames(top_genes)
write.csv(top_genes_list, "top_genes_list.csv")
saveRDS(vsd_corrected_ml, "vsd_corrected_ml.rds")
write.csv(vsd_corrected_ml, "vsd_corrected_ml.csv")


# Exploratory data analysis
# PCA
p_raw <- pca(vsd_raw, removeVar = 0.1)
p1 <- ggplot(p_raw$rotated, 
                     aes(x = PC1, y = PC2, color = factor(colData$batch))) + 
  geom_point(size = 1.5) + 
  labs(title = 'Before correction', color = "Batch") +
  theme_minimal()

p_corrected <- pca(vsd_corrected, removeVar = 0.1)
p2 <- ggplot(p_corrected$rotated,
                           aes(x = PC1, y = PC2, color = factor(colData$batch))) +
  geom_point(size = 1.5) +
  labs(title = 'After correction', color = "Batch") +
  theme_minimal()

# Print combined PCA plots
(p1 | p2) + plot_layout(guides = "collect") & theme(legend.position = "right")

# Volcano
resLFC_corrected <- lfcShrink(dds_corrected, coef = 'state_Th1_vs_Naive',
                              type = 'apeglm')
# Get gene symbols to switch with ENSEMBL IDs
resLFC_corrected$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(resLFC_corrected),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")
EnhancedVolcano(resLFC_corrected,
                lab = resLFC_corrected$symbol,
                title = "Volcano plot of DEGs",
                x = "log2FoldChange",
                y = "pvalue")

# Distance Heatmap
sample_dists <- dist(t(vsd_corrected))
sample_dist_matrix <- as.matrix(sample_dists)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         main = "Sample-to-sample distances")
                
                
                
                
                
                
                


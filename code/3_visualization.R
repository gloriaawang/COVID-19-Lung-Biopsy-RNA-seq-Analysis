# COVID-19 Lung Biopsy RNA-seq Analysis 
# Visualization Script

# Load required libraries
library(ggplot2)
library(gplots)
library(viridis)

# Read DESeq2 results
results <- read.csv("results/deseq2_results/all_results.csv", row.names=1)
norm_counts <- read.csv("results/deseq2_results/normalized_counts.csv", row.names=1)

# Create volcano plot
pdf("results/figures/volcano_plot.pdf", width=8, height=6)
plot(results$log2FoldChange, 
     -log10(results$padj),
     pch = 19,
     col = ifelse(results$padj < 0.05, "red", "black"),
     xlab = "log2 (Infected/Uninfected expression)",
     ylab = "-log10(adjusted p-value)",
     main = "")

# Add significance threshold line
abline(h = -log10(0.05), col = "green")

# Highlight IFI6 gene
IFI6_index <- which(rownames(results) %in% c('NM_022873', 'NM_002038'))
points(results$log2FoldChange[IFI6_index], 
       -log10(results$padj)[IFI6_index],
       pch = 19,
       col = "blue")
text(results$log2FoldChange[IFI6_index],
     -log10(results$padj)[IFI6_index],
     labels = "IFI6",
     pos = 4,
     col = "blue")
dev.off()

# Create expression heatmap
sig_genes <- subset(results, padj < 0.05)
sig_counts <- log2(norm_counts[rownames(sig_genes),] + 0.5)

pdf("results/figures/expression_heatmap.pdf", width=8, height=10)
heatmap.2(as.matrix(sig_counts),
          scale = "none",
          Colv = TRUE,
          Rowv = TRUE,
          col = viridis(50),
          trace = "none",
          density.info = "none",
          cexCol = 1,
          symkey = FALSE,
          symbreaks = FALSE,
          labRow = FALSE,
          main = "Gloria Wang",
          key.title = "",
          key.xlab = "log2 TPM")
dev.off()

# Create fold change heatmap
log2FC_matrix <- matrix(rep(sig_genes$log2FoldChange, 2), ncol = 2)
colnames(log2FC_matrix) <- c("log2FoldDiff", "log2FoldDiff2")

pdf("results/figures/fold_change_heatmap.pdf", width=6, height=10)
heatmap.2(log2FC_matrix,
          scale = "none",
          Colv = NA,
          Rowv = TRUE,
          col = viridis(50, option = "D"),
          trace = "none",
          density.info = "none",
          cexCol = 1,
          symkey = FALSE,
          symbreaks = FALSE,
          labRow = FALSE,
          main = "",
          key.title = "",
          key.xlab = "log2 Fold Change")
dev.off()

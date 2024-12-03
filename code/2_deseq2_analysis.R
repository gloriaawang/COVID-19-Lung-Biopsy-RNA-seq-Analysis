# COVID-19 Lung Biopsy RNA-seq Analysis
# Differential Expression Analysis

# Load required libraries
library(DESeq2)
library(dplyr)

# Read kallisto output files
read_abundance <- function(file) {
  read.table(file, header=TRUE, row.names=1)
}

# Import abundance files
uninfected1 <- read_abundance("kallisto_output/Uninfected_1/abundance.tsv")
uninfected2 <- read_abundance("kallisto_output/Uninfected_2/abundance.tsv")
infected <- read_abundance("kallisto_output/Infected/abundance.tsv")

# Create merged count matrix
counts <- data.frame(
  uninfected1 = uninfected1$est_counts,
  uninfected2 = uninfected2$est_counts,
  infected = infected$est_counts,
  row.names = rownames(uninfected1)
)

# Round counts to integers
counts <- round(counts)

# Create sample information
samples <- data.frame(
  'group' = c('uninfected', 'uninfected', 'infected'),
  row.names = c('uninfected1', 'uninfected2', 'infected')
)

# Create DESeq2 object and run analysis
dds <- DESeqDataSetFromMatrix(counts, samples, design = ~ group)
dds <- DESeq(dds)

# Get results for infected vs uninfected
infected_vs_uninfected <- results(dds, 
                                contrast = c('group', 'infected', 'uninfected'))

# Filter for significantly differentially expressed genes
significant_genes <- subset(infected_vs_uninfected, 
                          padj < 0.05)

# Save results
write.csv(as.data.frame(infected_vs_uninfected), 
          "results/deseq2_results/all_results.csv")
write.csv(as.data.frame(significant_genes),
          "results/deseq2_results/significant_genes.csv")

# Export normalized counts for visualization
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, 
          "results/deseq2_results/normalized_counts.csv")

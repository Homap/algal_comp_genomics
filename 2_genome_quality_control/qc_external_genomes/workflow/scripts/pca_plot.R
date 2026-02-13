library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript pca_plot.R <input.tsv> <output.pdf>")
}

input_file <- args[1]
output_file <- args[2]

df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Transform skewed variables
df <- df %>%
  mutate(
    log10_n50 = log10(n50 + 1),
    log10_num_contigs = log10(num_contigs + 1),
    log10_assembly_len = log10(assembly_len + 1)
  )

pca_data <- df %>%
  select(
    log10_n50,
    log10_num_contigs,
    log10_assembly_len,
    GC_content,
    gc_mean_per_contig,
    gc_sd_per_contig,
    N_content
  )

pca_res <- prcomp(pca_data, scale. = TRUE)

pca_scores <- as.data.frame(pca_res$x)
pca_scores$assembly <- df$assembly

# Variance explained (%)
var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
pc1_var <- round(var_explained[1], 1)
pc2_var <- round(var_explained[2], 1)

pdf(output_file)
ggplot(pca_scores, aes(x = PC1, y = PC2, label = assembly)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(
    title = "PCA of assembly statistics",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  )
dev.off()
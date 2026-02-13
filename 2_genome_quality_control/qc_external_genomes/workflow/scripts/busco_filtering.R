library(ggplot2)
library(dplyr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript pca_plot.R <input.tsv> <output.pdf>")
}

input_file <- args[1]
output_file <- args[2]

# load
busco <- fread(input_file)

# multivariate matrix
X <- as.matrix(busco[, .(C, C_D, F, M)])

set.seed(123)
km <- kmeans(X, centers = 2)

busco$cluster <- factor(km$cluster)
p <- prcomp(X)

pca_df <- as.data.frame(p$x)
pca_df$cluster <- factor(busco$cluster)
pca_df$sample <- busco$sample

# proportion of variance explained
var_explained <- (p$sdev^2) / sum(p$sdev^2)

pc1_perc <- round(var_explained[1] * 100, 1)
pc2_perc <- round(var_explained[2] * 100, 1)

ggplot(pca_df, aes(PC1, PC2, color = cluster)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = sample),
    size = 3,
    max.overlaps = Inf,   # show all labels
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
  theme_classic() +
  labs(
    x = paste0("PC1 (", pc1_perc, "%)"),
    y = paste0("PC2 (", pc2_perc, "%)")
  )


writeLines(cluster1_assemblies, output_file)

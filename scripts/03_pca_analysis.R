# ==============================================================================
# Script: 03_pca_analysis.R
# Description: Performs Principal Component Analysis (PCA) on the expression 
#              data and visualizes clustering by diet.
# Dependencies: Requires data/log_data.rds
# ==============================================================================
library(ggplot2)

log_data <- readRDS("data/log_data.rds")
numeric_data <- scale(t(log_data))

# ---------------------------------------------------------
# 1. FULL DATASET PCA
# ---------------------------------------------------------
pca_result <- prcomp(numeric_data, center = TRUE, scale. = TRUE)

num_samples <- ncol(log_data)
labels <- rep(1:(num_samples/3), each = 3)
shapes <- rep(c("triangle", "circle"), each = 3, times = num_samples / 6)
colors <- rep(1:(num_samples / 6), each = 6)

pca_data <- data.frame(
  Sample = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Diet = as.factor(labels),
  Shape = as.factor(shapes),
  ColorGroup = as.factor(colors)
)

# Plot and save Full PCA
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = ColorGroup, shape = Shape)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "PCA of Gene Expression by Diet (Full Dataset)",
       x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal() +
  scale_shape_manual(values = c(16, 17))

ggsave("outputs/03_pca_full_dataset.pdf", plot = pca_plot, width = 8, height = 6)

# ---------------------------------------------------------
# 2. EXTRACT SIGNIFICANT GENE LOADINGS (PC2)
# ---------------------------------------------------------
loadings <- pca_result$rotation[, 2]
mean_pc2 <- mean(loadings)
sd_pc2 <- sd(loadings)
n_sd <- 3 # Threshold of 3 standard deviations

significant_indices <- which(abs(loadings - mean_pc2) > n_sd * sd_pc2)

significant_genes_df <- data.frame(
  Gene = rownames(pca_result$rotation)[significant_indices],
  Loading = loadings[significant_indices]
)

# Order by absolute impact
significant_genes_df <- significant_genes_df[order(abs(significant_genes_df$Loading), decreasing = TRUE), ]

# Save top genes to a CSV for easy sharing/viewing
write.csv(significant_genes_df, "outputs/03_pca_significant_genes.csv", row.names = FALSE)
message("PCA analysis complete. Outputs saved to 'outputs/'.")
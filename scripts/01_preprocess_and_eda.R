# ==============================================================================
# Script: 01_preprocess_and_eda.R
# Description: Data loading, log-transformation, and initial exploratory data 
#              analysis (EDA) for primate diet gene expression.
# ==============================================================================

# Create directories if they don't exist
dir.create("outputs", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# ---------------------------------------------------------
# 1. LOAD AND TRANSFORM DATA
# ---------------------------------------------------------
url <- "https://raw.githubusercontent.com/ahmedmoustafa/gene-expression-datasets/refs/heads/main/datasets/primates_diet/primates_diet_clean.tsv"
raw_data <- read.table(url, header = TRUE, row.names = 1)
data_matrix <- as.matrix(raw_data)

# Log-transform the data
log_data <- log(data_matrix)

# ---------------------------------------------------------
# 2. EXPLORATORY DATA ANALYSIS (EDA)
# ---------------------------------------------------------
colors <- rep(c("pink", "red", "lavender", "purple", "lightgreen", "darkgreen", "lightblue", "navy"), each = 3)

# Save Boxplot
pdf("outputs/01_expression_boxplot.pdf", width = 10, height = 6)
boxplot(log_data, col = colors, main = "Log-Transformed Gene Expression", las = 2, cex.axis = 0.7)
dev.off()

# Save Initial Hierarchical Clustering
pdf("outputs/01_initial_clustering.pdf", width = 8, height = 6)
hc <- hclust(as.dist(1 - cor(log_data)))
plot(hc, main = "Hierarchical Clustering of Samples (Pre-filtering)")
dev.off()

# ---------------------------------------------------------
# 3. BATCH ORGANIZATION
# ---------------------------------------------------------
# Group columns into triplicates (samples)
matrix_list <- list()
for(i in seq(1, ncol(log_data), by = 3)) {
  max_index <- min(i + 2, ncol(log_data))
  matrix_list[[length(matrix_list) + 1]] <- log_data[, i:max_index]
}

# Separate into Batch 1 and Batch 2
matrix_list_batch1 <- matrix_list[seq(1, length(matrix_list), by = 2)]
matrix_list_batch2 <- matrix_list[seq(2, length(matrix_list), by = 2)]
combined_matrix_list <- c(matrix_list_batch1, matrix_list_batch2)

# Save the processed list for the next script
saveRDS(combined_matrix_list, "data/processed_matrices.rds")
saveRDS(log_data, "data/log_data.rds")
message("Preprocessing complete. Data saved to 'data/'.")
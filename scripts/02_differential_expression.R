# ==============================================================================
# Script: 02_differential_expression.R
# Description: Calculates fold changes and p-values, applies significance 
#              filters, and generates heatmaps.
# Dependencies: Requires data/processed_matrices.rds
# ==============================================================================
library(gplots)

matrix_list <- readRDS("data/processed_matrices.rds")
number_of_genes <- nrow(matrix_list[[1]])

# ---------------------------------------------------------
# 1. CALCULATE MEANS AND FOLD CHANGES
# ---------------------------------------------------------
means_list <- lapply(matrix_list, rowMeans)

folds_list <- list()
for(i in c(4, 8)) {
  for(j in (i-3):(i-1)) {
    folds_list[[length(folds_list) + 1]] <- list(i, j, means_list[[i]] - means_list[[j]])
  }
}

# ---------------------------------------------------------
# 2. STATISTICAL TESTING (T-TESTS)
# ---------------------------------------------------------
pvalues_list <- list()
for (j in 1:2) {
  for (k in 1:3) {
    pvalue <- numeric(number_of_genes)
    for (i in 1:number_of_genes) {
      x <- matrix_list[[4 * j]][i, ]
      y <- matrix_list[[k + 4 * (j - 1)]][i, ]
      
      if (sd(x) > 0 && sd(y) > 0) {
        pvalue[i] <- t.test(x, y)$p.value
      } else {
        pvalue[i] <- 0
      }
    }
    pvalues_list[[length(pvalues_list) + 1]] <- list(4 * j, k + 4 * (j - 1), pvalue)
  }
}

# ---------------------------------------------------------
# 3. FILTERING (VOLCANO PLOT LOGIC)
# ---------------------------------------------------------
fold_cutoff <- 1.0
pvalue_cutoff <- 0.05

pdf("outputs/02_volcano_plots.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))
for(i in 1:length(pvalues_list)) {
  batch <- ifelse(i <= 3, 1, 2)
  plot(folds_list[[i]][[3]], -log10(pvalues_list[[i]][[3]]),
       main = paste("Batch:", batch, "| Cat1 =", pvalues_list[[i]][[1]] - (batch - 1)*4, 
                    "Cat2 =", pvalues_list[[i]][[2]] - (batch - 1)*4),
       xlab = "Log2 Fold Change", ylab = "-Log10 P-value", pch = 16, col = "darkgray")
  
  abline(v = c(fold_cutoff, -fold_cutoff), col = "red", lwd = 2, lty = 2)
  abline(h = -log10(pvalue_cutoff), col = "blue", lwd = 2, lty = 2)
}
dev.off()

# Apply filters
filtered_list <- list()
for(i in 1:length(pvalues_list)) {
  sig_fold <- abs(folds_list[[i]][[3]]) >= fold_cutoff
  sig_pval <- pvalues_list[[i]][[3]] <= pvalue_cutoff
  combined_sig <- sig_fold & sig_pval
  
  j <- ifelse(i < 4, 0, 1)
  filtered_list[[i]] <- cbind(matrix_list[[i+j]][combined_sig, ], matrix_list[[4+4*j]][combined_sig, ])
}

# ---------------------------------------------------------
# 4. HEATMAP VISUALIZATION
# ---------------------------------------------------------
pdf("outputs/02_significant_heatmaps.pdf", width = 8, height = 8)
for(i in 1:length(filtered_list)) {
  batch <- ifelse(i <= 3, 1, 2)
  cat <- ifelse(i <= 3, i, i - 3)
  
  col_dend <- as.dendrogram(hclust(as.dist(1 - cor(filtered_list[[i]]))))
  row_dend <- as.dendrogram(hclust(as.dist(1 - cor(t(filtered_list[[i]])))))
  
  heatmap(filtered_list[[i]], Rowv = row_dend, Colv = col_dend, col = rev(redgreen(1024)),
          main = paste("Significant Genes - Batch:", batch, "| Cat:", cat))
}
dev.off()

saveRDS(filtered_list, "data/filtered_genes.rds")
message("Differential expression complete. Plots saved to 'outputs/'.")
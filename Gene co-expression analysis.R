# 2. Parameter optimization (finding optimal pK value)
# Scan pK values, pK parameter affects doublet detection sensitivity
sweep.res.list <- paramSweep(sc, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Visualize pK selection, choose pK value with highest BCmetric
ggplot(bcmvn, aes(x = as.factor(pK), y = BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "DoubletFinder Parameter Optimization", x = "pK", y = "BCmetric")

# Automatically select optimal pK value
optimal.pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
print(paste("Optimal pK value is:", optimal.pK))

# 3. Estimate doublet rate
# For 1200 cell dataset, doublet rate is typically 1% - 5%, here we use a common estimate
# You can also estimate more accurately based on experimental loading concentration
# Using correct doublet rate: about 0.8% for 1000 cells, about 1% for 1200 cells
nExp_poi <- round(0.01 * ncol(sc))  # 1% doublet rate, about 12 doublets for 1200 cells
print(paste("Using 1% doublet rate, estimated doublet count:", nExp_poi))

# Run DoubletFinder with optimal parameters pK=0.22
sc <- doubletFinder(sc, 
                    PCs = 1:30, 
                    pN = 0.25,
                    pK = 0.22,
                    nExp = nExp_poi, 
                    reuse.pANN = FALSE, 
                    sct = FALSE)

# View newly generated column names
new_cols <- grep("DF.classifications", names(sc@meta.data), value = TRUE)
doublet_col_name <- tail(new_cols, 1)
print(paste("Using doublet analysis column:", doublet_col_name))

# View overall results
doublet_info <- table(sc@meta.data[[doublet_col_name]])
print("Overall doublet classification results:")
print(doublet_info)

# Key analysis: Check doublet distribution in cluster 1 (290 cells)
cluster_cells <- WhichCells(sc, idents = 1)
cluster_doublet_info <- table(sc@meta.data[cluster_cells, doublet_col_name])

print(paste("Cluster 1 (total", length(cluster_cells), "cells) doublet distribution:"))
print(cluster_doublet_info)

# Calculate doublet count and proportion in cluster 1
if("Doublet" %in% names(cluster_doublet_info)) {
  doublet_count <- cluster_doublet_info["Doublet"]
  cluster_doublet_rate <- doublet_count / sum(cluster_doublet_info) * 100
  print(paste("Doublet count in cluster 1:", doublet_count, "cells"))
  print(paste("Doublet proportion in cluster 1:", round(cluster_doublet_rate, 2), "%"))
  
  # Key judgment
  if(doublet_count <= 2) {  # According to 0.8% rate, expect about 2-3 doublets in 290 cells
    print("✅ Doublet count in cluster 1 is within normal expected range")
    print("✅ This indicates cluster 1 is likely a real biological cell population")
  } else if(doublet_count > 10) {
    print("⚠️  Doublet count in cluster 1 is abnormally high")
    print("⚠️  This may indicate technical issues or need further inspection")
  } else {
    print("ℹ️  Doublet count in cluster 1 is within acceptable range")
  }
} else {
  print("✅ No doublets detected in cluster 1")
  print("✅ This strongly supports cluster 1 as a real homogeneous cell population")
}

# Visualization
p1 <- DimPlot(sc, reduction = "umap", group.by = doublet_col_name) +
  ggtitle(paste("Doublet Predictions (", nExp_poi, "expected doublets)")) +
  scale_color_manual(values = c("Singlet" = "grey", "Doublet" = "red"))
ggsave("Marker/doublet_cells.pdf", p1, width = 8, height = 6)
ggsave("Marker/doublet_cells.png", p1, width = 8, height = 6)

# Extract subset of cluster 1 cells
cluster1_cells <- WhichCells(sc, idents = 1)
sc_cluster1 <- subset(sc, cells = cluster1_cells)

# Draw scatter plot
p2 <- FeatureScatter(sc_cluster1, 
                     feature1 = "PTPRC", 
                     feature2 = "DCN",
                     pt.size = 1.5) +
  ggtitle("Single-cell level co-expression analysis: Cluster 1") +
  xlab("PTPRC expression level") +
  ylab("DCN expression level") +
  theme_minimal()
ggsave("Marker/DCN.pdf", p2, width = 8, height = 6)
ggsave("Marker/DCN.png", p2, width = 8, height = 6)
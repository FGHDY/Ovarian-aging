
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(batchelor)
library(SeuratWrappers)
library(magrittr)
library(tidyverse)
library(clusterProfiler)
library(GO.db)
library(org.Hs.eg.db)
library(DOSE)
library(DoubletFinder)
library(SingleR)
library(celldex)
library(harmony)
library(pheatmap)
library(openxlsx)
library(writexl)


setwd

scRNA<- readRDS("")

# Data normalization
scRNA <- NormalizeData(scRNA)  # Normalization
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)  # Select highly variable genes
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA))  # Standardization

str(scRNA)  # Check Seurat object structure

# View PCA results
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))  # PCA dimensionality reduction
plot4 <- DimPlot(scRNA, reduction = "pca",
                 group.by = "new.ident")  # Sample distribution in PCA
plot5 <- ElbowPlot(scRNA, ndims = 50, reduction = "pca")
# Draw elbow plot using first 30 dimensions

plot6 <- plot4 + plot5
dir.create("clustering")
ggsave("clustering/pca2.png", plot = plot6, width = 11, height = 6)
ggsave("clustering/pca2.pdf", plot = plot6, width = 11, height = 6)
pc.num = 1:30
# Cell clustering
scRNA <- FindNeighbors(scRNA, dims = pc.num)  # KNN+SNN algorithm
scRNA <- FindClusters(scRNA, resolution = 0.3)  # Louvain algorithm (optimal)
## clustree to select optimal resolution  
library(clustree)
scRNA <- FindClusters(scRNA, resolution = seq(0.1, 1.0, by = 0.1))
p1 <- clustree(scRNA@meta.data, prefix = "RNA_snn_res.")
colnames(scRNA@meta.data)
ggsave("clustering/bestresolution.png", plot = p1, width = 11, height = 7)

# Change cluster numbering from starting at 0 to starting at 1
scRNA@meta.data$seurat_clusters <- as.numeric(as.character(scRNA@meta.data$seurat_clusters)) + 1
scRNA@meta.data$seurat_clusters <- as.factor(scRNA@meta.data$seurat_clusters)

metadata <- scRNA@meta.data
# Save clustering results to CSV
cell_cluster <- data.frame(cell_ID = row.names(scRNA@meta.data),
                           cluster_ID = scRNA@meta.data$seurat_clusters)
write.csv(cell_cluster, 'clustering/cell_cluster.csv', row.names = F, quote = F)

# View cluster number distribution
table(scRNA@meta.data$seurat_clusters)

cluster_table <- table(scRNA@meta.data$seurat_clusters)  # Create table
write.csv(cluster_table, file = "ClusterTable.xlsx")  # Write table to xlsx

# Run TSNE
scRNA <- RunUMAP(scRNA, dims = pc.num)

# Draw UMAP plot (cluster numbering should now start from 1)
plot8 <- DimPlot(scRNA, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1)

# Display UMAP plot
print(plot7)
# View cluster distribution in UMAP dimensionality reduction
ggsave("clustering/UMAP1.png", plot = plot8, width = 9, height = 7)
ggsave("clustering/UMAP1.pdf", plot = plot8, width = 9, height = 7)

# View each sample's distribution in UMAP plot
plot8 = DimPlot(scRNA, reduction = "umap", group.by = 'new.ident', pt.size = 0.8)
# View cluster composition in each sample
plot9 = DimPlot(scRNA, reduction = "umap", split.by = 'new.ident',
                group.by = "seurat_clusters", ncol = 3, pt.size = 0.6, label = T)

ggsave("clustering/umap_new.ident.png", plot = plot8, width = 9, height = 7)
ggsave("clustering/umap_new.ident.pdf", plot = plot8, width = 9, height = 7)

ggsave("clustering/umap_new.ident2.png", plot9, width = 13, height = 7)
ggsave("clustering/umap_new.ident2.pdf", plot9, width = 13, height = 7)

saveRDS(scRNA, "umap_distribution.rds")









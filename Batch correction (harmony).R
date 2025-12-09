
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
scRNA <-readRDS("")

## Batch correction (harmony)
# harmony
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))  # PCA dimensionality reduction
scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "new.ident")
# Harmony requires PCA results; scRNA object already has PCA, not repeated
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dim = 1:30)  # UMAP dimensionality reduction

scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.2) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)  # Standard clustering

# Change cluster numbering from starting at 0 to starting at 1
scRNA_harmony@meta.data$seurat_clusters <- as.numeric(as.character(scRNA_harmony@meta.data$seurat_clusters)) + 1
scRNA_harmony@meta.data$seurat_clusters <- as.factor(scRNA_harmony@meta.data$seurat_clusters)

# Overall
p1 <- DimPlot(scRNA_harmony, group.by = "new.ident", pt.size = 0.1) +
  ggtitle("Integrated by harmony")
p2 <- DimPlot(scRNA, group.by = "new.ident", pt.size = 0.1) +
  ggtitle("No integrated")
p = p1 + p2 + plot_layout(guides = "collect")

ggsave("clustering/after_harmony.png", p1, width = 9, height = 7)
ggsave("clustering/after_harmony.pdf", p1, width = 9, height = 7)
# Each cluster
p3 <- DimPlot(scRNA_harmony, reduction = "umap", split.by = "new.ident", group.by = "seurat_clusters", label = T, ncol = 3) +
  ggtitle("Integrated by harmony")

p4 <- DimPlot(scRNA, reduction = "umap", split.by = "new.ident", group.by = "seurat_clusters", label = T, ncol = 3) +
  ggtitle("No integrated")
p5 <- p3 + p4 + plot_layout(guides = "collect")

ggsave("clustering/after_harmony3.png", p3, width = 13, height = 8)
ggsave("clustering/after_harmony3.png", p3, width = 13, height = 8)

p7 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "seurat_clusters", label = T) +
  ggtitle("Integrated by harmony")
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", label = T) +
  ggtitle("No integrated")
p9 <- p7 + p8 + plot_layout(guides = "collect")

ggsave("clustering/after_harmony2.png", p7, width = 9, height = 7)
ggsave("clustering/after_harmony2.pdf", p7, width = 9, height = 7)

## clustree to select optimal resolution  
library(clustree)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = seq(0.1, 1.0, by = 0.1))
clustree(scRNA_harmony@meta.data, prefix = "RNA_snn_res.")
colnames(scRNA@meta.data)
ggsave("clustering/bestresolution.png", p1, width = 11, height = 7)

saveRDS(scRNA_harmony, "harmony.rds")
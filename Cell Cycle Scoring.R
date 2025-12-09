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

setwd()
scRNA <- readRDS("")

#Cell Cycle Scoring

# 3.1 Load dataset - Human
cc_genes <- cc.genes.updated.2019  
s_genes <- intersect(cc_genes$s.genes, rownames(scRNA)) 
g2m_genes <- intersect(cc_genes$g2m.genes, rownames(scRNA))

## Mouse cell cycle genes
s_genes = c("Mcm5", "Pcna", "Tyms", "Mcm7", "Mcm4", "Rpa2", "Ung",
            "Clspn", "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2",
            "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2", "Cdc45")
g2m_genes = c("Hmgb2", "Cks2", "Tubb4b", "Hmgnb", "Cks1b", "Cdc20",
              "Tuba1b", "Top2a", "Ndc80", "Ckap2l", "Hmgb3", "Ckap5",
              "Aurkb", "Birc5", "Kif11", "Anp32e", "Tubb5", "Gtse1",
              "Kif20b", "Hmgb1", "Ckap2", "Cenpa", "Lbr", "Cks1brt",
              "Nuf2", "Ube2c", "Anln", "Cdk1", "Tacc3", "Fam64a", "Smc4")

# 3.2 Calculate cycle scores
scRNA_harmony <- CellCycleScoring(scRNA_harmony,
                                  s.features = s_genes,
                                  g2m.features = g2m_genes,
                                  set.ident = TRUE)

table(scRNA_harmony@meta.data$Phase)  # View phase distribution

# 3.3 Visualize cell cycle effect on PCA
scRNA_harmony <- RunPCA(scRNA_harmony, features = VariableFeatures(scRNA_harmony))
p1 <- DimPlot(scRNA, reduction = "pca", group.by = "Phase")  # Using scRNA object
ggsave("cell_cycle/pca.png", p1, width = 9, height = 5)

# Check PC correlation with cycle scores (<0.3 acceptable)
p2 <- FeatureScatter(scRNA, feature1 = "PC_1", feature2 = "S.Score", group.by = "Phase")
p3 <- FeatureScatter(scRNA, feature1 = "PC_2", feature2 = "G2M.Score", group.by = "Phase")
FeatureScatter(scRNA, feature1 = "S.Score", feature2 = "G2M.Score", group.by = "Phase")

ggsave("cell_cycle/cycle_correlation1.png", p2, width = 9, height = 5)
ggsave("cell_cycle/cycle_correlation2.png", p3, width = 9, height = 5)

## 4 Correct cell cycle effects
# 4.1 Regression to remove cycle effects
scRNA_harmony <- ScaleData(scRNA,
                           vars.to.regress = c("S.Score", "G2M.Score"),
                           features = rownames(scRNA))

# Re-run PCA for validation
scRNA_harmony <- RunPCA(scRNA_harmony, features = VariableFeatures(scRNA_harmony))
p1 <- DimPlot(scRNA_harmony, reduction = "pca", group.by = "Phase")
ggsave("cell_cycle/corrected_pca.png", p1, width = 9, height = 5)

## 5 Visualization and validation
# 5.1 UMAP comparison
# Before correction
scRNA_raw <- RunUMAP(scRNA_harmony, resolution = 0.2, dims = 1:30)  # Modified: scRNA
p1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "Phase", pt.size = 0.3) + ggtitle("Cell Cycle Phase")

ggsave("cell_cycle/before_correction_umap.png", p1, width = 9, height = 7)
ggsave("cell_cycle/before_correction_umap.pdf", p1, width = 9, height = 7)

# End of cell cycle analysis (no correction needed in normal studies)
saveRDS(scRNA_harmony, "cell_cycle.rds")

# After correction
scRNA_regressed <- RunUMAP(scRNA_harmony, resolution = 0.4, dims = 1:30)  # Note: using corrected scRNA
p2 <- DimPlot(scRNA_regressed, reduction = "umap", group.by = "Phase") + ggtitle("Cell Cycle Phase")

ggsave("cell_cycle/after_correction_umap.png", p2, width = 9, height = 5)
ggsave("cell_cycle/after_correction_umap.pdf", p2, width = 9, height = 5)

# 5.2 Differential gene validation
# Differential genes before correction
markers_raw <- FindMarkers(scRNA_raw, ident.1 = "Cluster1", ident.2 = "Cluster2")  # scRNA_raw

# Differential genes after correction
markers_regressed <- FindMarkers(scRNA_regressed, ident.1 = "Cluster1", ident.2 = "Cluster2")

## 6 Advanced analysis
# 6.1 Cyclone cycle prediction
library(scran)
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package = "scran"))

# Modified: using scRNA object
assignments <- cyclone(GetAssayData(scRNA, slot = "counts"), pairs = hs.pairs)
scRNA$Cyclone_Phase <- assignments$phases  # Store results in scRNA

saveRDS(scRNA_regressed, "cell_cycle_corrected.rds")
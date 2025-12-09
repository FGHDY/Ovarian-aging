### Cell Trajectory Analysis
#### CYTOTRACE2 Analysis
## Load R packages
library(CytoTRACE2)
library(Seurat)
library(tidyverse)
library(paletteer)
library(BiocParallel)

# Use Seurat object
sc <- readRDS("cell_subtype_classification.rds")

# Data preprocessing
cytotrace2_res <- cytotrace2(sc,  # Seurat object
                             is_seurat = TRUE,
                             slot_type = "counts",  # Both counts and data are acceptable
                             species = 'human')  # Species selection, default is mouse

# Visualization
annotation <- data.frame(phenotype = sc@meta.data$celltype) %>%
  set_rownames(., colnames(sc))

## plotData generates multiple plots at once, then stores in a list, use $ to view
plots <- plotData(cytotrace2_result = cytotrace2_res,
                  annotation = annotation,
                  is_seurat = TRUE)

p1 <- plots$CytoTRACE2_UMAP + theme(axis.text = element_text(size = 15))

p2 <- plots$Phenotype_UMAP +
  ggtitle("Celltype") +
  theme(
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15) )

p4 <- plots$CytoTRACE2_Relative_UMAP + theme(axis.text = element_text(size = 15))
p5 <- plots$CytoTRACE2_Boxplot_byPheno + theme(axis.text = element_text(size = 13)) +
  ggtitle("Phenotype") + theme(plot.title = element_text(size = 13))

p6 <- p1 + p2 + p4 + p5 + 
  plot_layout(
    ncol = 2, 
    guides = "collect",  # Unify collection and management of legends
    heights = c(1, 1)    # Uniform row height
  ) +
  theme(
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 15),
    legend.box = "vertical",  # Vertical arrangement of legends
    legend.spacing.y = unit(0.1, "cm")  # Adjust legend spacing
  )

ggsave("monocle/cytotrace1.png", p5, width = 5, height = 3.5)
ggsave("monocle/cytotrace1.pdf", p5, width = 5, height = 3.5)
saveRDS(cytotrace2_res, "cytotrace.rds")

### Monocle
# Load R packages
library(monocle)
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(ggsci)
library(igraph)
library(data.table)
library(Biobase)
library(Matrix)

sc <- readRDS("reclustering.rds")
### 2 Construct cds object by extracting data from Seurat object as an example

## Load Seurat information to construct Monocle object
expr_matrix <- as(as.matrix(sc@assays$RNA@layers$counts), 'sparseMatrix')

expr_matrix <- sc@assays$RNA@layers$counts  # Use sparse matrix directly to avoid conversion to dense matrix
# Extract table information to p_data (phenotype_data)
p_data <- sc@meta.data
p_data$celltype <- sc@active.ident  ## Integrate cell identification information for each cell into p_data
## Extract gene information, such as biotype, GC content, etc.
f_data <- data.frame(gene_short_name = row.names(sc), row.names = row.names(sc))
## The number of rows in expr_matrix is the same as the number of rows in f_data (gene number), and the number of columns in expr_matrix is the same as the number of columns in p_data.

## Construct CDS object CellDataSet
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
# Convert p_data and f_data from data.frame to AnnotatedDataFrame objects

cds <- newCellDataSet(
  as(expr_matrix, "sparseMatrix"),
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size())

## 3 Estimate size.factor and dispersion
cds <- estimateSizeFactors(cds, locfunc = median)
cds <- estimateDispersions(cds)
cds <- estimateDispersions(
  cds,
  method = "pooled",          # Speed up calculation
  fitType = "local",          # Avoid complex fitting
  min_cells_detected = 10,    # Only process genes expressed in ≥10 cells
  cores = 4                   # Multi-core parallel (if available)
)

## 4 Filter low-quality cells (optional, as Seurat QC has been done)
cds <- detectGenes(cds, min_expr = 0.1)  # This will add a column num_cells_expressed to fdata
# View cell phenotype data
head(pData(cds))

# View gene annotation data
head(fData(cds))

expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))
cds <- cds[expressed_genes, ]
# The above step filters out genes expressed in less than 10 cells.

## 5 Cell ordering
# Monocle website provides 4 selection methods: 1. Select developmentally differentially expressed genes; 2. Select cluster differentially expressed genes;
# 3. Select genes with high dispersion; 4. Custom developmental marker genes
# The first three are unsupervised analysis methods, and cell developmental trajectory generation is completely without manual intervention;
# The last is a semi-supervised analysis method that can use prior knowledge
# 5.1 Select genes that define the process
seurat1 = sc
## Use Seurat to select highly variable genes (not recommended)
express_genes <- VariableFeatures(seurat1)
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)

## Use cluster differentially expressed genes (recommended)
# Stricter differential gene screening
deg.clusters <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# express_genes <- subset(deg.clusters, p_val_adj < 0.05)$gene  # Not used, directly use method 1 below

# A better method is to select genes that vary between different cell states (recommended!)
# Method 1: Select genes that are differentially expressed in multiple clusters
top_genes_per_cluster <- deg.clusters %>% group_by(cluster) %>% top_n(200, avg_log2FC)
selected_genes <- unique(top_genes_per_cluster$gene)

## Use Monocle to select highly variable genes
# Method 2: Directly use genes with high dispersion (more suitable for Monocle)
disp_genes <- dispersionTable(cds)
disp_genes <- disp_genes %>% arrange(desc(dispersion_empirical))
selected_genes <- head(disp_genes$gene_id, 2000)

# Set ordering genes
cds <- setOrderingFilter(cds, selected_genes)
plot_ordering_genes(cds)

# 5.2 Dimensionality reduction: Monocle2 uses DDRTree algorithm
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

# 5.3 Pseudotime axis trajectory construction and cell ordering in pseudotime
# Order cells and complete trajectory construction based on the expression trend of ordering genes
cds <- orderCells(cds)

# Use root_state parameter to set the root of the pseudotime axis. As can be seen from the pseudotime coloring diagram below, the left side is the root.
# According to the state diagram, the root is State1. If you want to set the other end as the root, you can do the following:
# cds <- orderCells(cds, root_state = 5)  # Set State5 as the starting point of the pseudotime axis

# 5.4 Visualization: Color cells based on phenotype information (metadata) in cds@phenoData@data
# 1) Color by pseudotime value (Pseudotime is the probability calculated by Monocle2 based on cell gene expression information, representing the sequence in time.)
dir.create("monocle")

# If celltype names in the data are incorrect
# Ensure cell names match (important!)
rownames(pData(cds)) <- colnames(cds)  # Ensure row names are consistent

# Extract correct celltype information from original sc object
pData(cds)$celltype <- sc$celltype[match(rownames(pData(cds)), colnames(sc))]

# Verify fix result
table(pData(cds)$celltype)  # Should now show GC1-GC5

p1 <- plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
ggsave("monocle/Pseudotime.png", p1, width = 9, height = 6)
ggsave("monocle/Pseudotime.pdf", p1, width = 9, height = 6)
# cds <- orderCells(cds, root_state = 5)  # Set State5 as the starting point of the pseudotime axis
dev.off()  # Close the current graphics device and reset the graphics device list to the initial state.

cds <- orderCells(cds, root_state = 5)
# 2) Color by cell type
p2 <- plot_cell_trajectory(cds, color_by = "celltype",
                           size = 1, label.size = 2, show_backbone = TRUE) +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 5,           # Dot size
        alpha = 1,          # Transparency
        stroke = 0.5         # Border thickness
      ),
      nrow = 1,             # Number of legend columns
      keywidth = unit(1, "cm"),   # Legend key width
      keyheight = unit(0.6, "cm") # Legend key height
    )
  )

ggsave("monocle/celltype.png", p2, width = 9, height = 6)
ggsave("monocle/celltype.pdf", p2, width = 9, height = 6)
dev.off()

# 3) Color by cell state
p3 <- plot_cell_trajectory(cds, color_by = "State", size = 3, label = TRUE, show_backbone = TRUE) +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 5,           # Dot size
        alpha = 1,          # Transparency
        stroke = 0.5         # Border thickness
      ),
      nrow = 1,             # Number of legend columns
      keywidth = unit(1, "cm"),   # Legend key width
      keyheight = unit(0.6, "cm") # Legend key height
    )
  )

ggsave("monocle/State.png", p3, width = 9, height = 6)
ggsave("monocle/State.pdf", p3, width = 9, height = 6)
dev.off()

# 4) Order cells by Seurat clusters (consistent with previous celltype classification results)
p4 <- plot_cell_trajectory(cds, color_by = "seurat_clusters", size = 3, label = TRUE, show_backbone = TRUE)

dev.off()

# 5) Color by cell state (split) trajectory diagram to more easily see the position of each state in 3).
p5 <- plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 2) +
  theme(text = element_text(size = 16))
ggsave("monocle/State_split.png", p5, width = 17, height = 10)                                                              
dev.off()
# Combine plots
p <- p1 + p2 + p3
ggsave("monocle/combined_plot.png", p, width = 17, height = 10)   

# 4. orig.ident
p4 <- plot_cell_trajectory(cds, color_by = 'new.ident') +
  scale_color_manual(values = allcolour) + 
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 5,           # Dot size
        alpha = 1,          # Transparency
        stroke = 0.5         # Border thickness
      ),
      nrow = 1,             # Number of legend columns
      keywidth = unit(1, "cm"),   # Legend key width
      keyheight = unit(0.6, "cm") # Legend key height
    )
  )
p4
ggsave("monocle/new.ident.png", p4, width = 9, height = 6)
ggsave("monocle/new.ident.pdf", p4, width = 9, height = 6)

# Pseudotime heatmap
diff_test_res <- differentialGeneTest(cds,
                                      reducedModelFormulaStr = "~new.ident",
                                      cores = 4)

gene_to_cluster <- diff_test_res %>% 
  arrange(qval) %>% 
  head(50) %>% 
  pull(gene_short_name)
head(gene_to_cluster)

p1 <- plot_pseudotime_heatmap(cds[gene_to_cluster, ],
                              num_clusters = 5, 
                              show_rownames = TRUE,
                              cores = 4, return_heatmap = TRUE,
                              hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))

ggsave("monocle/pseudotime_heatmap1.png", p1, width = 7, height = 9)
ggsave("monocle/pseudotime_heatmap1.pdf", p1, width = 7, height = 9)

# Gene trajectory plot
gene <- head(gene_to_cluster)
gene <- c("RPS5", "RPS2", 'RPS3')
p2 <- plot_cell_trajectory(cds, markers = gene,
                           use_color_gradient = TRUE)

ggsave("monocle/HSP.png", p2, width = 17, height = 10)
ggsave("monocle/HSP.pdf", p2, width = 17, height = 10)

# Gene pseudotime scatter plot
gene <- c("HSPB1", "CXCL12", "RUNDC3B")

plot_genes_in_pseudotime(cds[gene, ],
                         color_by = "new.ident",
                         nrow = 3,
                         ncol = NULL)
ggsave("monocle/gene_scatter_plot2.png", p3, width = 7, height = 7)
ggsave("monocle/gene_scatter_plot2.pdf", p3, width = 7, height = 7)

saveRDS(cds, "monocle2_analysis.rds")

cds <- readRDS("monocle2_analysis.rds")
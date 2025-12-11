#Data Consolidation Analysis
# Load required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(cowplot)
library(stringr)

setwd("")

# 1. Read conserved gene table
conserved_genes <- read.csv("D:/ovaryaging/三者联合分析/三者同源基因/conserved_genes.csv", 
                            stringsAsFactors = FALSE)

# View data structure
convert_gene_names_filtered <- function(seurat_obj, gene_map, species_name) {
  cat("Processing", species_name, "dataset (filtering empty strings)...\n")
  
  # Get expression matrix
  expr_data <- GetAssayData(seurat_obj, slot = "counts")
  cat("Original gene count:", nrow(expr_data), "\n")
  
  # Find genes in mapping table
  current_genes <- rownames(expr_data)
  genes_in_map <- intersect(current_genes, names(gene_map))
  cat("Genes in mapping table:", length(genes_in_map), "\n")
  
  if (length(genes_in_map) == 0) {
    stop("Error: No genes found in mapping table!")
  }
  
  # Convert to new gene names
  new_gene_names <- gene_map[genes_in_map]
  
  # Filter out empty gene names
  non_empty_idx <- !is.na(new_gene_names) & new_gene_names != ""
  genes_in_map <- genes_in_map[non_empty_idx]
  new_gene_names <- new_gene_names[non_empty_idx]
  
  cat("Gene count after filtering empty strings:", length(new_gene_names), "\n")
  
  if (length(new_gene_names) == 0) {
    stop("Error: All gene names are empty strings!")
  }
  
  # Create mapping data frame
  mapping_df <- data.frame(
    original = genes_in_map,
    new = new_gene_names,
    stringsAsFactors = FALSE
  )
  
  # Check for duplicates
  if (any(duplicated(new_gene_names))) {
    cat("Found duplicate genes, removing duplicates...\n")
    
    # Identify duplicate genes
    dup_genes <- new_gene_names[duplicated(new_gene_names)]
    cat("Number of duplicate genes:", length(unique(dup_genes)), "\n")
    
    # Show first 5 duplicate genes
    if (length(dup_genes) > 0) {
      cat("First 5 duplicate genes:\n")
      print(head(unique(dup_genes), 5))
    }
    
    # Keep only first occurrence of each gene
    keep_idx <- !duplicated(new_gene_names)
    mapping_df <- mapping_df[keep_idx, ]
  }
  
  # Create new expression matrix
  new_expr <- expr_data[mapping_df$original, , drop = FALSE]
  rownames(new_expr) <- mapping_df$new
  
  cat("Final gene count:", nrow(new_expr), "\n")
  
  # Create new Seurat object
  new_obj <- CreateSeuratObject(
    counts = new_expr,
    project = species_name,
    meta.data = seurat_obj@meta.data
  )
  
  return(new_obj)
}

# First check for empty values in conserved gene table
cat("Checking for empty values in conserved gene table...\n")

# Check for empty strings in each column
empty_human <- sum(conserved_genes$Human.gene.name == "" | is.na(conserved_genes$Human.gene.name))
empty_mouse <- sum(conserved_genes$Mouse.gene.name == "" | is.na(conserved_genes$Mouse.gene.name))
empty_goat <- sum(conserved_genes$Goat.gene.name == "" | is.na(conserved_genes$Goat.gene.name))

cat("Empty values in human gene names:", empty_human, "\n")
cat("Empty values in mouse gene names:", empty_mouse, "\n")
cat("Empty values in goat gene names:", empty_goat, "\n")

# View rows with empty values
if (empty_mouse > 0) {
  cat("\nRows with empty mouse gene names:\n")
  empty_rows <- which(conserved_genes$Mouse.gene.name == "" | is.na(conserved_genes$Mouse.gene.name))
  print(conserved_genes[empty_rows, ])
  
  # Remove rows with empty gene names
  conserved_genes_clean <- conserved_genes[!conserved_genes$Mouse.gene.name == "" & !is.na(conserved_genes$Mouse.gene.name), ]
  cat("\nRows in cleaned conserved gene table:", nrow(conserved_genes_clean), "\n")
  
  # Recreate gene mappings
  mouse_gene_map_clean <- setNames(conserved_genes_clean$Human.gene.name, conserved_genes_clean$Mouse.gene.name)
  human_gene_map_clean <- setNames(conserved_genes_clean$Human.gene.name, conserved_genes_clean$Human.gene.name)
  goat_gene_map_clean <- setNames(conserved_genes_clean$Human.gene.name, conserved_genes_clean$Goat.gene.name)
  
  # Use cleaned mappings
  mouse_gene_map <- mouse_gene_map_clean
  human_gene_map <- human_gene_map_clean
  goat_gene_map <- goat_gene_map_clean
}

# Use filtered function
mouse_converted <- convert_gene_names_filtered(mouse, mouse_gene_map, "Mouse")
human_converted <- convert_gene_names_filtered(human, human_gene_map, "Human")
goat_converted <- convert_gene_names_filtered(goat, goat_gene_map, "Goat")

# Continue with integration analysis
cat("\n=== Continuing integration analysis ===\n")

# 1. Check converted data
cat("Checking dimensions of converted data:\n")
cat("Human:", dim(human_converted), "\n")
cat("Mouse:", dim(mouse_converted), "\n")
cat("Goat:", dim(goat_converted), "\n")

# 2. Find common genes
cat("\nFinding common genes across three species...\n")
common_genes <- Reduce(intersect, list(
  rownames(human_converted),
  rownames(mouse_converted),
  rownames(goat_converted)
))
cat("Number of common genes across three species:", length(common_genes), "\n")

# 3. Keep only common genes
cat("Keeping only common genes...\n")
human_converted <- subset(human_converted, features = common_genes)
mouse_converted <- subset(mouse_converted, features = common_genes)
goat_converted <- subset(goat_converted, features = common_genes)

cat("Dimensions after processing:\n")
cat("Human:", dim(human_converted), "\n")
cat("Mouse:", dim(mouse_converted), "\n")
cat("Goat:", dim(goat_converted), "\n")

# 4. Add metadata
cat("\nAdding metadata...\n")

# Human
human_converted$species <- "Human"
human_converted$age <- human_converted$new.ident
human_converted$celltype <- human_converted$celltype
human_converted$group <- paste("Human", human_converted$age, sep = "_")
human_converted$species_age <- paste("Human", human_converted$age, sep = "_")

# Mouse
mouse_converted$species <- "Mouse"
mouse_converted$age <- mouse_converted$new.ident
mouse_converted$celltype <- mouse_converted$celltype
mouse_converted$group <- paste("Mouse", mouse_converted$age, sep = "_")
mouse_converted$species_age <- paste("Mouse", mouse_converted$age, sep = "_")

# Goat
goat_converted$species <- "Goat"
goat_converted$age <- goat_converted$new.ident
goat_converted$celltype <- goat_converted$celltype
goat_converted$group <- paste("Goat", goat_converted$age, sep = "_")
goat_converted$species_age <- paste("Goat", goat_converted$age, sep = "_")

# 5. Create object list
cat("\nCreating object list...\n")
obj.list <- list(human_converted, mouse_converted, goat_converted)
names(obj.list) <- c("Human", "Mouse", "Goat")

# 6. Data normalization and feature selection
cat("\nData normalization and feature selection...\n")
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(obj.list[[i]], 
                                 normalization.method = "LogNormalize",
                                 scale.factor = 10000,
                                 verbose = FALSE)
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]],
                                        selection.method = "vst",
                                        nfeatures = 3000,
                                        verbose = FALSE)
}

# 7. Select integration features
cat("\nSelecting integration features...\n")
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
cat("Number of selected integration features:", length(features), "\n")

# Revised integration workflow
cat("\n=== Revised integration workflow ===\n")

# 1. Re-run normalization and feature selection
cat("\n1. Re-running normalization and feature selection...\n")
for (i in 1:length(obj.list)) {
  cat("Processing", names(obj.list)[i], "...\n")
  
  # Ensure using RNA assay
  DefaultAssay(obj.list[[i]]) <- "RNA"
  
  # Check if normalization already done
  if (!"RNA" %in% Assays(obj.list[[i]]) || 
      is.null(obj.list[[i]][["RNA"]]@data) || 
      all(dim(obj.list[[i]][["RNA"]]@data) == 0)) {
    cat("  Running NormalizeData...\n")
    obj.list[[i]] <- NormalizeData(
      obj.list[[i]],
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  } else {
    cat("  Normalized data already exists\n")
  }
  
  # Check if variable features already found
  var_features <- VariableFeatures(obj.list[[i]])
  if (is.null(var_features) || length(var_features) == 0) {
    cat("  Running FindVariableFeatures...\n")
    obj.list[[i]] <- FindVariableFeatures(
      obj.list[[i]],
      selection.method = "vst",
      nfeatures = 2000,
      verbose = FALSE
    )
  } else {
    cat("  Variable features already exist, count:", length(var_features), "\n")
  }
}

# 2. Select integration features
cat("\n2. Selecting integration features...\n")
features <- SelectIntegrationFeatures(
  object.list = obj.list,
  nfeatures = 2000
)
cat("Number of selected features:", length(features), "\n")

# 3. Run ScaleData and PCA
cat("\n3. Running ScaleData and PCA...\n")
for (i in 1:length(obj.list)) {
  cat("Processing", names(obj.list)[i], "...\n")
  
  # Run ScaleData
  cat("  Running ScaleData...\n")
  obj.list[[i]] <- ScaleData(
    obj.list[[i]],
    features = features,
    verbose = FALSE
  )
  
  # Run PCA
  cat("  Running PCA...\n")
  obj.list[[i]] <- RunPCA(
    obj.list[[i]],
    features = features,
    npcs = 30,
    verbose = FALSE
  )
  
  # Check PCA results
  if (!"pca" %in% Reductions(obj.list[[i]])) {
    cat("  Warning: PCA not successful, trying with all genes...\n")
    obj.list[[i]] <- RunPCA(
      obj.list[[i]],
      npcs = 30,
      verbose = FALSE
    )
  }
  
  cat("  PCA dimensions:", dim(Embeddings(obj.list[[i]], "pca")), "\n")
}

# 4. Find integration anchors
cat("\n4. Finding integration anchors...\n")
tryCatch({
  anchors <- FindIntegrationAnchors(
    object.list = obj.list,
    anchor.features = features,
    reduction = "rpca",
    dims = 1:30,
    k.anchor = 20,
    verbose = TRUE
  )
  cat("Anchor finding successful!\n")
}, error = function(e) {
  cat("rpca failed, trying cca...\n")
  anchors <- FindIntegrationAnchors(
    object.list = obj.list,
    anchor.features = features,
    reduction = "cca",  # Fall back to cca method
    dims = 1:30,
    k.anchor = 20,
    verbose = TRUE
  )
})

# 9. Integrate data
cat("\nIntegrating data...\n")
combined <- IntegrateData(
  anchorset = anchors,
  dims = 1:30,
  verbose = TRUE
)

# 10. Downstream analysis
cat("\nPerforming downstream analysis...\n")

# Set default assay
DefaultAssay(combined) <- "integrated"

# Normalization
combined <- ScaleData(combined, verbose = FALSE)

# PCA
combined <- RunPCA(combined, 
                   npcs = 50, 
                   verbose = FALSE,
                   features = features)

# Check PCA results
cat("\nChecking PCA results...\n")
pca_plot <- DimPlot(combined, 
                    reduction = "pca", 
                    group.by = "species",
                    pt.size = 0.1) +
  ggtitle("PCA by Species")
print(pca_plot)

# Elbow plot
elbow_plot <- ElbowPlot(combined, ndims = 50) +
  ggtitle("PCA Elbow Plot")
print(elbow_plot)

# Select number of PCs (based on elbow plot)
n.pcs <- 30
cat("Using", n.pcs, "PCs for dimensionality reduction\n")

# 11. Run UMAP
cat("\nRunning UMAP...\n")
combined <- RunUMAP(
  combined,
  reduction = "pca",
  dims = 1:n.pcs,
  n.neighbors = 30,
  min.dist = 0.3,
  metric = "cosine",
  verbose = FALSE
)

# 12. Run t-SNE (optional)
cat("\nRunning t-SNE (optional)...\n")
combined <- RunTSNE(combined,
                    dims = 1:n.pcs,
                    perplexity = 30,
                    verbose = FALSE)

# 13. Clustering
cat("\nPerforming cell clustering...\n")
combined <- FindNeighbors(combined, dims = 1:n.pcs)
combined <- FindClusters(combined, resolution = c(0.1, 0.2, 0.5, 1.0))

# Use appropriate resolution
Idents(combined) <- "integrated_snn_res.0.5"
cat("Clustering results (resolution 0.5):\n")
print(table(Idents(combined)))

# 14. Define color schemes
cat("\nDefining color schemes...\n")

# Species colors
species_colors <- c("Human" = "#E41A1C", "Mouse" = "#377EB8", "Goat" = "#4DAF4A")

# Age colors
age_colors <- c("young" = "#66C2A5", "old" = "#FC8D62")

# Group colors
group_colors <- c("Human_young" = "#E41A1C", "Human_old" = "#F781BF",
                  "Mouse_young" = "#377EB8", "Mouse_old" = "#A6CEE3",
                  "Goat_young" = "#4DAF4A", "Goat_old" = "#B2DF8A")

# Cell type colors (based on your provided cell types)
celltype_colors <- c(
  "SC" = "#1F77B4",    # Granulosa cells (possibly)
  "GC" = "#FF7F0E",    # Granulosa cells
  "TC" = "#2CA02C",    # Theca cells
  "LC" = "#D62728",    # Oocytes
  "EC" = "#9467BD",    # Endothelial cells
  "SMC" = "#8C564B",   # Smooth muscle cells
  "IC" = "#E377C2",    # Immune cells
  "T" = "#7F7F7F",     # T cells
  "Ma" = "#BCBD22",    # Macrophages
  "B" = "#17BECF",     # B cells
  "LEC" = "#9EDAE5",   # Lymphatic endothelial cells
  "EpC" = "#FF9896",   # Epithelial cells
  "SwC" = "#C5B0D5"    # Stromal cells
)

# 15. Visualization function
create_umap_plot <- function(object, group.by, colors = NULL, title = "", 
                             pt.size = 0.1, label = FALSE, label.size = 3) {
  p <- DimPlot(object, 
               reduction = "umap", 
               group.by = group.by,
               cols = colors,
               pt.size = pt.size,
               label = label,
               label.size = label.size,
               repel = TRUE,
               raster = FALSE) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          legend.text = element_text(size = 8),
          legend.title = element_blank())
  return(p)
}

# 16. Create basic visualizations
cat("\nCreating basic visualizations...\n")

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
# By species
p1 <- create_umap_plot(combined, "species", species_colors, "By Species")+
  ggtitle("") +
  theme(legend.text = ggtext::element_markdown(size = 18),  # Using package prefix
        legend.title = element_text(size = 16))
dir.create("cluster")
ggsave("cluster/uamp1.png",p1,width = 9,height = 6)
ggsave("cluster/uamp1.pdf",p1,width = 9,height = 6)

# By age
p2 <- create_umap_plot(combined, "age", age_colors, "By Age") + ggtitle("") +
  theme(legend.text = ggtext::element_markdown(size = 18),  # Using package prefix
        legend.title = element_text(size = 16))
ggsave("cluster/uamp2.png",p2,width = 9,height = 6)
ggsave("cluster/uamp2.pdf",p2,width = 9,height = 6)

# By species and age combination
p3 <- create_umap_plot(combined, "group", group_colors, "By Species and Age")+
  ggtitle("") +
  theme(legend.text = ggtext::element_markdown(size = 18),  # Using package prefix
        legend.title = element_text(size = 16))

ggsave("cluster/uamp3.png",p3,width = 9,height = 6)
ggsave("cluster/uamp3.pdf",p3,width = 9,height = 6)

# By cell type
p4 <- create_umap_plot(combined, "celltype", celltype_colors, "By Cell Type", label.size = 8,label = TRUE)+
  ggtitle("") +
  theme(legend.text = ggtext::element_markdown(size = 18),  # Using package prefix
        legend.title = element_text(size = 16))
ggsave("cluster/uamp4.png",p4,width = 9,height = 6)
ggsave("cluster/uamp4.pdf",p4,width = 9,height = 6)
# By cluster
p5 <- create_umap_plot(combined, "integrated_snn_res.0.5", NULL, "By Cluster (res=0.5)", label = TRUE)

# Combined plot
library(patchwork)
combined_plot <- (p1 + p2) / (p3 + p4) / (p5 + plot_spacer()) + 
  plot_annotation(tag_levels = 'A')
print(combined_plot)

# 17. Faceted visualizations
cat("\nCreating faceted visualizations...\n")

# Adjust species order, age order
# Convert species column to factor and set level order
combined$species <- factor(combined$species, 
                           levels = c("Human", "Mouse", "Goat"))
combined$age <- factor(combined$age,
                       levels = c("young","old"))
combined$species_age <- factor(combined$species_age,
                               levels = c("Human_young","Mouse_young","Goat_young",
                                          "Human_old","Mouse_old","Goat_old"))
combined$celltype <- factor(combined$celltype,
                            levels = c("SC","GC","LC","EC","T","TC","SMC","IC","Ma","EpC","LEC","B","SwC"))


# Facet by species
p_facet_species <- DimPlot(combined,
                           reduction = "umap",
                           split.by = "species",
                           group.by = "celltype",
                           ncol = 3,
                           label = T,
                           label.size = 8,
                           pt.size = 0.5,
                           cols = celltype_colors) +
  ggtitle("Cell Types by Species") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_facet_species)
ggsave("cluster/uamp5.png",p_facet_species,width = 16,height = 6)
ggsave("cluster/uamp5.pdf",p_facet_species,width = 16,height = 6)
# Facet by age
p_facet_age <- DimPlot(combined,
                       reduction = "umap",
                       split.by = "age",
                       group.by = "celltype",
                       ncol = 2,
                       pt.size = 0.1,
                       cols = celltype_colors) +
  ggtitle("Cell Types by Age") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_facet_age)
ggsave("cluster/uamp6.png",p_facet_age,width = 10,height = 6)
ggsave("cluster/uamp6.pdf",p_facet_age,width = 10,height = 6)


# Facet by species and age combination
p_facet_combined <- DimPlot(combined,
                            reduction = "umap",
                            split.by = "species_age",
                            group.by = "celltype",
                            ncol = 3,
                            pt.size = 1,
                            cols = celltype_colors) +
  ggtitle("Cell Types by Species and Age") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

print(p_facet_combined)
ggsave("cluster/uamp7.png",p_facet_combined,width = 16,height = 8)
ggsave("cluster/uamp7.pdf",p_facet_combined,width = 16,height = 8)
# 18. Cell type composition analysis
cat("\nPerforming cell type composition analysis...\n")

# Create statistics table
library(dplyr)
cell_stats <- combined@meta.data %>%
  group_by(species, age, celltype) %>%
  summarize(n_cells = n(), .groups = 'drop') %>%
  group_by(species, age) %>%
  mutate(percent = round(n_cells / sum(n_cells) * 100, 2)) %>%
  ungroup()

# Print statistics
cat("\nCell type statistics:\n")
print(cell_stats)

# 19. Visualize cell type composition
library(ggplot2)

# Add combined label text
cell_stats <- cell_stats %>%
  mutate(label_text = paste0(n_cells, "\n(", percent, "%)"))

p_bar <- ggplot(cell_stats, aes(x = celltype, y = percent, fill = age)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ species, ncol = 3) +
  scale_fill_manual(values = age_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "top",
        strip.text = element_text(face = "bold")) +
  labs(x = "Cell Type", y = "Percentage (%)",
       title = "Cell Type Composition by Species and Age",
       fill = "Age") +
  geom_text(aes(label = label_text),  # Use combined label
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 2.5)

print(p_bar)

# 20. Heatmap: Cell type proportions
cat("\nCreating cell type proportion heatmap...\n")
library(pheatmap)

# Prepare heatmap data
cell_wide <- cell_stats %>%
  select(species, age, celltype, percent) %>%
  pivot_wider(names_from = c(species, age), 
              values_from = percent, 
              values_fill = 0)

# Convert to matrix
cell_matrix <- as.matrix(cell_wide[, -1])
rownames(cell_matrix) <- cell_wide$celltype

# Create heatmap
P <- pheatmap(cell_matrix,
              main = "Cell Type Percentage Heatmap",
              color = colorRampPalette(c("white", "red"))(100),
              display_numbers = TRUE,
              number_format = "%.1f",
              cluster_rows = TRUE,
              cluster_cols = FALSE,
              fontsize_row = 19,
              fontsize_col = 12,
              fontsize_number = 12,
              angle_col = 45)
ggsave("cluster/heatmap1.png",P,width = 10,height = 6)
ggsave("cluster/heatmap1.pdf",P,width = 10,height = 6)
# 21. Marker gene analysis
cat("\nPerforming marker gene analysis...\n")

# Set default assay to RNA for differential expression analysis
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined, normalization.method = "LogNormalize")

# Find cell type marker genes
Idents(combined) <- "celltype"
# After normalization and scaling, use JoinLayers function to ensure data layers are correctly associated with metadata layers.
combined<- JoinLayers(combined, prefix = "RNA")
celltype_markers <- FindAllMarkers(combined,
                                   only.pos = TRUE,
                                   min.pct = 0.25,
                                   logfc.threshold = 0.25,
                                   test.use = "wilcox")

# Save marker gene results
write.csv(celltype_markers, "celltype_markers.csv", row.names = FALSE)
cat("Cell type marker genes saved to: celltype_markers.csv\n")

# Get top 5 marker genes for each cell type
top_markers <- celltype_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Create heatmap
DefaultAssay(combined) <- "RNA"
combined <- ScaleData(combined, features = rownames(combined))

heatmap_plot <- DoHeatmap(combined,
                          features = top_markers$gene,
                          group.by = "celltype",
                          group.colors = celltype_colors,
                          size = 3,
                          angle = 0) +
  ggtitle("Top Marker Genes for Each Cell Type") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(heatmap_plot)
dir.create("Marker")
ggsave("Marker/heatmap.png",heatmap_plot,width = 22,height = 15)
ggsave("Marker/heatmap.pdf",heatmap_plot,width = 22,height = 15)
# 22. Cross-species comparison analysis
cat("\nPerforming cross-species comparison analysis...\n")

# Compare expression of same cell type across different species
# For example, compare granulosa cells (GC) across species
DefaultAssay(combined) <- "RNA"

# Create subset of GC cells
gc_cells <- subset(combined, subset = celltype == "GC")

# Compare differences between species
if (ncol(gc_cells) > 0) {
  Idents(gc_cells) <- "species"
  gc_markers <- FindAllMarkers(gc_cells,
                               only.pos = TRUE,
                               min.pct = 0.1,
                               logfc.threshold = 0.25)
  
  cat("Number of marker genes for GC cells across species:", nrow(gc_markers), "\n")
  
  # Save results
  write.csv(gc_markers, "gc_celltype_across_species_markers.csv", row.names = FALSE)
}

# 23. Age-related change analysis
cat("\nPerforming age-related change analysis...\n")

# Analyze age-related differential expression for each species
age_deg_results <- list()

for (sp in c("Human", "Mouse", "Goat")) {
  cat("Analyzing age-related differential expression for", sp, "...\n")
  
  # Extract data for this species
  sp_cells <- subset(combined, subset = species == sp)
  
  if (ncol(sp_cells) > 10) {  # Ensure sufficient cells
    Idents(sp_cells) <- "age"
    
    # Find differentially expressed genes between young vs old
    age_markers <- FindMarkers(sp_cells,
                               ident.1 = "old",
                               ident.2 = "young",
                               min.pct = 0.1,
                               logfc.threshold = 0.25)
    
    # Add gene name column
    age_markers$gene <- rownames(age_markers)
    age_markers$species <- sp
    
    age_deg_results[[sp]] <- age_markers
    
    # Save results
    write.csv(age_markers, paste0(sp, "_age_DEGs.csv"), row.names = FALSE)
  }
}

# Combine age-related DEGs for all species
all_age_degs <- do.call(rbind, age_deg_results)
write.csv(all_age_degs, "all_species_age_DEGs.csv", row.names = FALSE)

# 24. Save results
cat("\nSaving integration results...\n")

# Save integrated Seurat object
saveRDS(combined, file = "integrated_ovary_3species_final.rds")
cat("Integrated object saved to: integrated_ovary_3species_final.rds\n")

# Save UMAP coordinates
umap_coords <- as.data.frame(Embeddings(combined, reduction = "umap"))
tsne_coords <- as.data.frame(Embeddings(combined, reduction = "tsne"))
metadata <- combined@meta.data

all_data <- cbind(metadata, umap_coords)
write.csv(all_data, 
          file = "integrated_umap_coordinates_final.csv",
          row.names = TRUE)
cat("UMAP coordinates saved to: integrated_umap_coordinates_final.csv\n")

# 25. Statistical summary
cat("\n=== Integration Analysis Statistical Summary ===\n")
cat("Total cell count:", ncol(combined), "\n")
cat("Total gene count:", nrow(combined), "\n")
cat("\nStatistics by species:\n")
print(table(combined$species))
cat("\nStatistics by age:\n")
print(table(combined$age))
cat("\nStatistics by species and age:\n")
print(table(combined$species, combined$age))
cat("\nCell type distribution:\n")
print(table(combined$celltype))
cat("\nCell type distribution by species:\n")
print(table(combined$celltype, combined$species))

# 26. Create HTML report (optional)
cat("\nCreating HTML report (optional)...\n")
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  # Create simple report
  report_text <- c(
    "---",
    "title: 'Cross-Species Ovarian Single-Cell Integration Analysis Report'",
    "output: html_document",
    "---",
    "",
    "## Analysis Overview",
    paste0("- Analysis date: ", Sys.Date()),
    paste0("- Total cell count: ", ncol(combined)),
    paste0("- Total gene count: ", nrow(combined)),
    paste0("- Number of species: 3 (Human, Mouse, Goat)"),
    "",
    "## Cell Statistics",
    "```{r echo=FALSE}",
    "library(knitr)",
    "kable(table(combined$species, combined$age), caption = 'Cell Count Statistics')",
    "kable(table(combined$celltype), caption = 'Cell Type Distribution')",
    "```",
    "",
    "## UMAP Visualization",
    "### By Species",
    "```{r echo=FALSE, fig.width=10, fig.height=6}",
    "DimPlot(combined, reduction = 'umap', group.by = 'species', pt.size = 0.1)",
    "```",
    "",
    "### By Cell Type",
    "```{r echo=FALSE, fig.width=12, fig.height=8}",
    "DimPlot(combined, reduction = 'umap', group.by = 'celltype', pt.size = 0.1, label = TRUE)",
    "```"
  )
  
  writeLines(report_text, "integration_report.Rmd")
  # rmarkdown::render("integration_report.Rmd")  # Uncomment to generate report
}

# 27. Save all images
cat("\nSaving all images...\n")

# Save combined plot
ggsave("combined_umap_plots.png", combined_plot, 
       width = 16, height = 20, dpi = 300, bg = "white")

# Save faceted plots
ggsave("facet_by_species.png", p_facet_species, 
       width = 18, height = 6, dpi = 300, bg = "white")

ggsave("facet_by_age.png", p_facet_age, 
       width = 12, height = 6, dpi = 300, bg = "white")

ggsave("celltype_composition.png", p_bar, 
       width = 16, height = 8, dpi = 300, bg = "white")

# Save heatmap
png("celltype_markers_heatmap.png", width = 1200, height = 800, res = 150)
print(heatmap_plot)
dev.off()

# 28. Create interactive visualizations (optional)
cat("\nCreating interactive visualizations (optional)...\n")
if (requireNamespace("plotly", quietly = TRUE)) {
  library(plotly)
  
  # 3D UMAP
  combined <- RunUMAP(combined, dims = 1:n.pcs, n.components = 3, reduction.name = "umap3d")
  umap3d_coords <- Embeddings(combined, "umap3d")
  
  # Create interactive 3D plot
  plot_3d <- plot_ly(x = umap3d_coords[,1], 
                     y = umap3d_coords[,2], 
                     z = umap3d_coords[,3],
                     type = "scatter3d", 
                     mode = "markers",
                     color = combined$species,
                     colors = species_colors,
                     marker = list(size = 2),
                     text = paste("Cell type:", combined$celltype, 
                                  "<br>Species:", combined$species,
                                  "<br>Age:", combined$age),
                     hoverinfo = "text")
  
  # Save as HTML
  htmlwidgets::saveWidget(plot_3d, file = "interactive_3d_umap.html")
  cat("Interactive 3D UMAP saved to: interactive_3d_umap.html\n")
}

# 29. Session information
cat("\n=== Session Information ===\n")
sessionInfo()

cat("\n=== Analysis Complete! ===\n")
cat("Main output files:\n")
cat("1. Integrated object: integrated_ovary_3species_final.rds\n")
cat("2. UMAP coordinates: integrated_umap_coordinates_final.csv\n")
cat("3. Cell type marker genes: celltype_markers.csv\n")
cat("4. Age-related differential genes: all_species_age_DEGs.csv\n")
cat("5. Visualization images: various PNG format images\n")
cat("\nNext step suggestions:\n")
cat("1. Check cell type separation in UMAP plots\n")
cat("2. Validate biological significance of marker genes\n")
cat("3. Perform deeper cell subpopulation analysis\n")
cat("4. Perform cell communication analysis\n")

saveRDS(combined,"combined.rds")
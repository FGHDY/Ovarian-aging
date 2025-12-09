
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
scRNA_harmony <-readRDS("")

Marker gene identification
# Identify marker genes
dir.create("./Marker")

# Switch dataset idents from clusters back to orig.ident (samples)
Idents(scRNA_harmony) <- 'seurat_clusters'
# After normalization and standardization, use JoinLayers to ensure data layers are correctly associated with metadata layers
scRNA_harmony <- JoinLayers(scRNA_harmony, prefix = "RNA")

# Identify marker genes for each cluster, with selected cluster as experimental group and all remaining clusters as control group
all.markers = FindAllMarkers(scRNA_harmony, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25, 
                             only.pos = TRUE)
head(all.markers)

# Sort by FoldChange to select top10 marker genes for each cell cluster
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
# Replace n with any number; if only want top1, replace 10 with 1
top3 = all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

library(openxlsx)
write.xlsx(all.markers, 
           "Marker/all_Markers_of_each_clusters.xlsx", 
           colnames = T, 
           rownames = F, 
           sep = "\t")
write.xlsx(top10, 
           "Marker/top10_Markers_of_each_clusters.xlsx", 
           colnames = T, 
           rownames = F, 
           sep = "\t")
write.xlsx(top20, 
           "Marker/top20_Markers_of_each_clusters.xlsx", 
           colnames = T, 
           rownames = F, 
           sep = "\t")
write.xlsx(top50, 
           "Marker/top50_Markers_of_each_clusters.xlsx", 
           colnames = T, 
           rownames = F, 
           sep = "\t")

# Draw heatmap
scRNA_harmony <- ScaleData(scRNA_harmony, features = row.names(scRNA_harmony))  # Standardize all genes as differential genes may appear outside highly variable genes
heatmap_plot = DoHeatmap(object = scRNA_harmony, 
                         features = as.vector(top10$gene),  # Visualize top10 markers for each cluster
                         group.by = "seurat_clusters",  # Group by column
                         group.bar = T,  # Draw colorbar indicating cluster positions
                         size = 2) +
  theme(axis.text.y = element_text(size = 4))  # Adjust text size

# Save plot results
ggsave("Marker/top10_marker_of_each_cluster_heatmap.png", width = 24, height = 12,
       plot = heatmap_plot)
ggsave("Marker/top10_marker_of_each_cluster_heatmap.pdf", width = 24, height = 12,
       plot = heatmap_plot)
heatmap_plot 



# Draw violin plot for single gene
VlnPlot(scRNA_harmony, 
        cols = my36colors, 
        features = c("GENE1"),
        pt.size = 0.1,  # Choose not to display points representing cells; can change to non-zero value (integer or decimal) to see cell expression level distribution
        group.by = "seurat_clusters")  # Display by cluster


## Grouped by new.ident
VlnPlot(scRNA_harmony, 
        cols = my36colors, 
        features = c("GENE2"),
        pt.size = 0.1,  # Choose not to display points representing cells; can change to non-zero value (integer or decimal) to see cell expression level distribution
        group.by = "new.ident")  # Display by cluster


### Cell annotation
### Apply annotation results to object's meta.data for subsequent advanced analysis
dir.create("singleR")

### Marker reference dataset bubble plot
markers <- unique(c(
  "DCN", "PDGFRA", "LUM",
  "DCN", "PDGFRA", "LUM",
  "HSD3B2", "CYP19A1", "HSD11B1",
  "DCN", "PDGFRA", "LUM",
  "DCN", "PDGFRA", "LUM",
  "VWF", "FLT1", "PECAM1",
  "HSD3B2", "CYP19A1", "HSD11B1",
  "RGS6", "ITGA7", "NOTCH3",
  "CD163", "PTPRC", "CD53",
  "INHBB", "AMH", "FSHR",
  "MMRN1", "FLT4", "CD36",
  "INHBB", "AMH", "FSHR",
  "KRT19", "PAX8", "CDH1",
  "KRT19", "PAX8", "CDH1",
  "STAR", "CYP17A1", "CYP11A1",
  "CDH19", "SOX10", "LGI4"
))

scRNA_harmony <- RenameIdents(scRNA_harmony,
                              "1" = "SC",
                              "2" = "SC",
                              "3" = "LC",
                              "4" = "SC",
                              "5" = "SC",
                              "6" = "EC",
                              "7" = "LC",
                              "8" = "SMC",
                              "9" = "IC",
                              "10" = "GC",
                              "11" = "LEC",
                              "12" = "GC",
                              "13" = "EpC",
                              "14" = "EpC",
                              "15" = "TC",
                              "16" = "SwC")

scRNA_harmony[["celltype"]] <- Idents(scRNA_harmony)
head(scRNA_harmony@meta.data$seurat_clusters)
color <- c('#6CB8D2', '#3C77AF', '#7DBFA7', '#EE934E', '#D1352C',
           '#B383B9', '#BBDD78', '#8FA4AE', '#479D88',
           '#415284', '#E69F84', '#D0AFC4', '#95538B', '#46976D')

table(scRNA_harmony$celltype)

# Visualization
p1 <- DotPlot(scRNA_harmony, features = markers, group.by = "celltype",
              dot.scale = 8,
              cols = c("blue", "red")) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.title.x = element_blank(),  # Remove X-axis title
        axis.title.y = element_blank(),  # Remove Y-axis title
        axis.text.y = element_text(size = 18)) +  # Adjust axis text size
  labs(x = NULL, y = NULL)  # Ensure axis titles are empty

ggsave("SingleR/cell_marker_genes.png", p1, width = 11, height = 6, bg = "white")
ggsave("SingleR/cell_marker_genes.pdf", p1, width = 11, height = 6, bg = "white")


#### Cell proportion histogram
table(scRNA_harmony$celltype)  # View cell counts per group
prop.table(table(Idents(scRNA_harmony)))
table(Idents(scRNA_harmony), scRNA_harmony$new.ident)  # Cell counts per group and cell cluster
## Generate cross-group cell count table
count_table <- table(Idents(scRNA_harmony), scRNA_harmony$new.ident)
# Save as CSV file
write.csv(count_table, "SingleR/cell_counts_by_group.csv")

Cellratio <- prop.table(table(Idents(scRNA_harmony), scRNA_harmony$new.ident), margin = 2)  # Calculate proportions of different cell clusters per sample group
Cellratio <- as.data.frame(Cellratio)  ## Convert format; if using percentages later, don't convert yet

# Save as CSV file
write.csv(Cellratio, "SingleR/cell_ratio_by_group.csv")

# Convert to percentage values (multiply by 100)
Cellratio_percent <- Cellratio * 100
Cellratio_percent <- as.data.frame(Cellratio_percent)
# Save as CSV file
write.csv(Cellratio_percent, "SingleR/cell_ratio_by_group.csv")

allcolour = c("#DC143C", "#0000FF", "#20B2AA", "#FFA500", "#9370DB", "#98FB98", "#F08080", "#1E90FF", "#7CFC00", "#FFFF00",
              "#808000", "#FF00FF", "#FA8072", "#7B68EE", "#9400D3", "#800080", "#A0522D", "#D2B48C", "#D2691E", "#87CEEB", "#40E0D0", "#5F9EA0",
              "#FF1493", "#0000CD", "#008B8B", "#FFE4B5", "#8A2BE2", "#228B22", "#E9967A", "#4682B4", "#32CD32", "#F0E68C", "#FFFFE0", "#EE82EE",
              "#FF6347", "#6A5ACD","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
library(ggplot2)
             
# Assume Cellratio is the correct proportion table and has been calculated by row (sample groups)
p5 <- ggplot(Cellratio_percent, aes(x = Var2, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity", width = 0.7, linewidth = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = 'Sample', y = 'Percentage(%)') +
  scale_fill_manual(values = color) + theme(axis.title.x = element_blank(),
                                            axis.text = element_text(size = 15),
                                            axis.title.y = element_text(size = 14),
                                            legend.text = element_text(size = 16),
                                            legend.title = element_blank(),
                                            strip.text = element_text(size = 20),
                                            legend.spacing.y = unit(0.5, "cm"),
                                            legend.key.size = unit(1,"cm"))

ggsave("SingleR/proportion_annotation2.png", p5, width = 4, height = 6)
ggsave("SingleR/proportion_annotation2.pdf", p5, width = 4, height = 6)


# Calculate raw cell counts for different cell clusters per group
Cellcount <- table(Idents(scRNA_harmony), scRNA_harmony$new.ident)
Cellcount <- as.data.frame(Cellcount)
library(ggplot2)

# Draw bar plot using raw cell counts
p6 <- ggplot(Cellcount, aes(x = Var2, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity", width = 0.7, linewidth = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = 'Sample', y = 'Cell Count') +
  scale_fill_manual(values = color) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        strip.text = element_text(size = 20),
        legend.spacing.y = unit(0.5, "cm"),
        legend.key.size = unit(1,"cm")) +
  geom_text(aes(label = sprintf("%d", Freq)), 
            position = position_stack(vjust = 0.5), 
            vjust = 0.5, 
            color = 'black')

ggsave("SingleR/count_annotation.png", p6, width = 5, height = 7)
ggsave("SingleR/count_annotation.pdf", p6, width = 5, height = 7)

# Group display
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales) # For percentage formatting
library(ggrepel)

# 1. Count cell type proportions by group
grouped_counts <- scRNA_harmony@meta.data %>%
  group_by(new.ident, celltype) %>%  # Group by both group and cell type
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)  # Percentage within group

# 2. Grouped histogram
P3 <- ggplot(grouped_counts, 
             aes(x = celltype, y = percentage, fill = new.ident)) +  # Fill colored by group
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  # Display side by side
  scale_fill_manual(values = c("young" = "#2c7bb6", "old" = "#d7191c")) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 3) +
  scale_y_continuous(labels = percent_format(scale = 1)) +  # Y-axis shows percentage
  labs(x = "Cell Type", y = "Percentage", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Faceted pie chart
p4 <- ggplot(grouped_counts, aes(x = "", y = percentage, fill = celltype)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  geom_text(aes(label = paste0(round(percentage, 0.5), "%")), 
            position = position_stack(vjust = 0.5)) +
  facet_wrap(~ new.ident) +  # Facet by group
  labs(fill = "Cell Type") +
  theme_void()

ggsave("SingleR/cell_proportion2.pdf", P3, width = 12, height = 5)
ggsave("SingleR/cell_proportion2.png", P3, width = 12, height = 5)

ggsave("SingleR/cell_proportion.pdf", p4, width = 12, height = 6)
ggsave("SingleR/cell_proportion.png", p4, width = 12, height = 6)

library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(scales)
library(grid)
# Create proportion data frame
prop_data <- data.frame(
  CellType = c("SC", "LC", "EC", "SMC", "IC", "GC", "LEC", "EpC", "TC", "SwC"),
  young = c(0.4818753878, 0.2974386245, 0.0708588141, 0.0215368253, 
            0.0296906851, 0.0654524506, 0.0191438447, 0.0033678986, 
            0.0098821235, 0.0007533457),
  old = c(0.8370957239, 0.0272411548, 0.0350553506, 0.0359778598,
          0.0246906881, 0.0026047319, 0.0117212937, 0.0234968526,
          0.0003255915, 0.0017907532))
# Convert to long format
df_long <- prop_data %>%
  pivot_longer(cols = -CellType, 
               names_to = "Group", 
               values_to = "Proportion") %>%
  mutate(Proportion = Proportion * 100,  # Convert to percentage
         Group = factor(Group, levels = c("young", "old")))

# Create independent charts for each cell type
create_cell_plot <- function(cell_type) {
  # Filter data for specific cell type
  cell_df <- df_long %>% filter(CellType == cell_type)
  
  # Calculate change and determine color
  change <- (cell_df$Proportion[cell_df$Group == "old"] - 
               cell_df$Proportion[cell_df$Group == "young"])
  
  line_color <- ifelse(change > 0, "#d7191c", "#2c7bb6")  # Red for increase, blue for decrease
  
  # Create chart
  p <- ggplot(cell_df, aes(x = Group, y = Proportion, group = CellType)) +
    geom_line(color = line_color, size = 1.5) +
    geom_point(aes(fill = Group), size = 5, shape = 21, color = "white", stroke = 1.5) +
    geom_text(aes(label = sprintf("%.1f%%", Proportion)), 
              vjust = -1.5, size = 4, fontface = "bold") +
    scale_fill_manual(values = c("young" = "#2c7bb6", "old" = "#d7191c")) +
    labs(title = cell_type,
         y = "Percentage (%)",
         x = NULL) +
    ylim(0, max(cell_df$Proportion) * 1.3) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_blank(),
      axis.title = element_text(size = 12),
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#f9f9f9", color = NA),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.8),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # Add change arrow and label
  if (change != 0) {
    arrow_x <- ifelse(change > 0, 1.5, 1.5)
    p <- p + 
      annotate("segment", 
               x = 1, xend = 2, 
               y = cell_df$Proportion[1], yend = cell_df$Proportion[2],
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = line_color, size = 1) +
      annotate("text", x = 1.5, y = mean(cell_df$Proportion),
               label = paste0(ifelse(change > 0, "+", ""), sprintf("%.1f%%", abs(change))),
               color = line_color, size = 4, fontface = "bold")
  }
  
  return(p)
}

# Create charts for each cell type
cell_plots <- lapply(prop_data$CellType, create_cell_plot)

# Remove all Y-axis titles except the first plot
for (i in 1:length(cell_plots)) {
  cell_plots[[i]] <- cell_plots[[i]] + 
    theme(axis.title.y = element_blank())
}

# Arrange all charts in a grid
combined_plot <- grid.arrange(
  grobs = cell_plots,
  ncol = 4,
  left = textGrob("Percentage(%)", 
                  gp = gpar(fontsize = 14, fontface = "bold"),
                  rot = 90),  # Y-axis title displayed vertically
  top = textGrob("human", 
                 gp = gpar(fontsize = 16, fontface = "bold")),
  padding = unit(1, "cm")
)

ggsave("SingleR/cell_type_proportion_change.png", combined_plot,
       width = 12, height = 7.5, dpi = 300, bg = "white")
ggsave("SingleR/cell_type_proportion_change.pdf", combined_plot,
       width = 12, height = 7.5)

saveRDS(scRNA_harmony, "cell_classification.rds")


## Enrichment Analysis
# Identify marker genes for each cluster, with selected cluster as experimental group and all remaining clusters as control group

all.markers = FindAllMarkers(sc, 
                             min.pct = 0.2, 
                             logfc.threshold = 0.25, 
                             only.pos = TRUE)  # TRUE extracts only positively expressed genes
head(all.markers)
# Replace n with any number; if only want top1, replace 10 with 1
top5 = all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
library(org.Hs.eg.db)
library(clusterProfiler)
top100 = all.markers
# Convert SYMBOL to ENTREZID
ids = bitr(top100$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db')
top100 = merge(top100, ids, by.x = 'gene', by.y = 'SYMBOL')

## KEGG annotation
gcSample = split(top100$ENTREZID, top100$cluster)
# KEGG
KEGG <- compareCluster(gcSample,
                       fun = "enrichKEGG",
                       organism = "hsa",
                       pvalueCutoff = 0.05)
p1 <- dotplot(KEGG)

p1

ggsave("enrichment_analysis_old/KEGG.pdf", p1, width = 5, height = 7)
ggsave("enrichment_analysis_old/KEGG.png", p1, width = 5, height = 7)

# GO enrichment
GO <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05)

p2 <- dotplot(GO) 

p2

ggsave("enrichment_analysis_old/GO.pdf", p2, width = 5, height = 7)
ggsave("enrichment_analysis_old/GO.png", p2, width = 5, height = 7)


















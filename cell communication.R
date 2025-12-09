

## Load R packages
library(CellChat)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(NMF)
library(ggalluvial)
library(patchwork)

# Set working directory
setwd("D:/ovaryaging/three_species_analysis/human_aging_ovary_new/GSE202601_RAW")

# Load data
scRNA <- readRDS("cell_classification.rds")

# Set ident object
Idents(scRNA) <- 'new.ident'

### Construct CellChat object
## Extract data needed for CellChat object and construct manually (not recommended)
data.input <- GetAssayData(scRNA, assay = "RNA", slot = "data")  # Normalized data matrix
labels <- data.frame(celltype = labels, row.names = as.character(labels$celltype))
# Create a dataframe
# Construct CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")

### Directly construct CellChat object from Seurat object (automatic method, recommended)
cellchat <- createCellChat(object = scRNA, group.by = "celltype")
levels(cellchat@idents)  # View cell clusters

# Assume scRNA is your Seurat object
# Extract young group cells
young_cells <- subset(scRNA, subset = new.ident == "young")
# Extract old group cells
old_cells <- subset(scRNA, subset = new.ident == "old")

# Create CellChat object for young group
cellchat_young <- createCellChat(young_cells, 
                                 group.by = "celltype")

# Create CellChat object for old group
cellchat_old <- createCellChat(old_cells, 
                               group.by = "celltype")

### Set database
# Set ligand-receptor interaction database
CellChatDB <- CellChatDB.human  # CellChatDB.mouse/CellChatDB.zebrafish
# Use CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB
# Use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")  # Secreted signaling
# ECM-Receptor (Extracellular matrix-receptor)
# Cell-Cell Contact (Cell-cell contact interactions)
#### Add database to CellChat object ####
cellchat@DB <- CellChatDB.use

cellchat_young@DB <- CellChatDB.use  # Young group
cellchat_old@DB <- CellChatDB.use     # Old group

# View database composition using human as an example
showDatabaseCategory(CellChatDB.human)

### Data preprocessing
# Extract subset of expression data cellchat@data.signaling
cellchat <- CellChat::subsetData(cellchat)

cellchat_young <- CellChat::subsetData(cellchat_young)  # Young group
cellchat_old <- CellChat::subsetData(cellchat_old)     # Old group

# View expression matrix
cellchat@data.signaling[1:10, 1:10]

cellchat_young@data.signaling[1:10, 1:10]  # Young
cellchat_old@data.signaling[1:10, 1:10]    # Old

## Identify overexpressed genes cellchat@var.features
cellchat <- identifyOverExpressedGenes(cellchat)

cellchat_young <- identifyOverExpressedGenes(cellchat_young)
cellchat_old <- identifyOverExpressedGenes(cellchat_old)

# View data structure
str(cellchat@var.features)

#### Identify overexpressed ligand-receptor pairs cellchat@LR ####
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat_young <- identifyOverExpressedInteractions(cellchat_young)
cellchat_old <- identifyOverExpressedInteractions(cellchat_old)

## View data structure
str(cellchat@LR)

#### Calculate cell communication probability
# Calculate communication probability based on ligand-receptor pairs, and calculate p-value based on permutation test,
# Infer communication network, and store in cellchat@net.
cellchat <- computeCommunProb(cellchat)

cellchat_young <- computeCommunProb(cellchat_young)
cellchat_old <- computeCommunProb(cellchat_old)

# Filter out communications between cell groups with too few cells
cellchat <- filterCommunication(cellchat, min.cells = 50)

cellchat_young <- filterCommunication(cellchat_young, min.cells = 10)
cellchat_old <- filterCommunication(cellchat_old, min.cells = 10)

## Infer cell-cell communication at pathway level and store in cellchat@netP
cellchat <- computeCommunProbPathway(cellchat)

cellchat_young <- computeCommunProbPathway(cellchat_young)
cellchat_old <- computeCommunProbPathway(cellchat_old)

## Extract all communication networks as dataframe format
# Can change slot.name
df.net <- subsetCommunication(cellchat, slot.name = "netP")

df.net.young <- subsetCommunication(cellchat_young, slot.name = "netP")
df.net.old <- subsetCommunication(cellchat_old, slot.name = "netP")

head(df.net)

## Calculate aggregated cell communication network ####
#'object@net$count' is a matrix: rows and columns are sources and targets respectively, 
#'and elements are the number of interactions between any two cell groups. 
#'USER can convert a matrix to a data frame using the function 'reshape2::melt()'
#'object@net$weight' is also a matrix containing the interaction weights between any two cell groups
cellchat <- aggregateNet(cellchat)

cellchat_young <- aggregateNet(cellchat_young)
cellchat_old <- aggregateNet(cellchat_old)

# Save CellChat objects
saveRDS(cellchat, "cellchat.rds")

saveRDS(cellchat_young, "cellchat_young.rds")
saveRDS(cellchat_old, "cellchat_old.rds")

### Visualize cell communication

# Use circle plot to show number of interactions or interaction strength between two cell groups
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd = TRUE)

netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Interaction weights/strength")
table(cellchat@meta$new.ident)

# Young
groupSize <- as.numeric(table(cellchat_young@idents))
par(mfrow = c(1, 2), xpd = TRUE)
# Create first plot (interaction count)
netVisual_circle(cellchat_young@net$count,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Number of interactions") 

# Create second plot (interaction weight)
netVisual_circle(cellchat_young@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Interaction weights/strength")  

table(cellchat_young@meta$new.ident)

dir.create("cellchat")

# Old
groupSize <- as.numeric(table(cellchat_old@idents))
par(mfrow = c(1, 2), xpd = TRUE)

netVisual_circle(cellchat_old@net$count,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Number of interactions")
netVisual_circle(cellchat_old@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Interaction weights/strength")

table(cellchat_old@meta$new.ident)

### Circle plot separate display
mat <- cellchat@net$weight
par(mfrow = c(2, 3), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = TRUE, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}

cellchat@netP$pathways
pathways.show <- c("EGF") 
levels(cellchat@meta$celltype)
# Hierarchy plot
vertex.receiver = seq(3, 4)
netVisual_aggregate(cellchat, signaling = pathways.show, 
                    vertex.receiver = vertex.receiver, 
                    layout = "hierarchy")

# Network plot Circle plot: netVisual_circle
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
# par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", weight.size = 6,
                    height = 8)

# Heatmap
par(mfrow = c(1, 1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

### Cell communication network
# Calculate network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# Heatmap visualization
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show,
                                  width = 12, height = 6)
# Signaling role analysis on aggregated cell communication network
netAnalysis_signalingRole_scatter(cellchat)
# Heatmap
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")

ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

p6 = ht1 + ht2

### Communication patterns
library(NMF)
library(ggalluvial)
# Identify key signals and potential communication patterns in all signaling pathways, k.range sets x-axis range
selectK(cellchat, pattern = "incoming", k.range = seq(2, 8))  # Or "incoming"
# Based on the number of patterns in the above figure, choose pattern <= 4 when number is highest

nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming",
                                          k = nPatterns,
                                          height = 10)
# River plot
netAnalysis_river(cellchat, pattern = "incoming")

# Dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

# Violin plot
plotGeneExpression(cellchat, signaling = "EGF", 
                   type = "violin")  # signaling can be changed arbitrarily
# Bubble plot
plotGeneExpression(cellchat, signaling = "EGF", 
                   type = "dot", 
                   color.use = "red")  # signaling can be changed arbitrarily

### Customization
# Bubble plot (overall)
# The x-axis of the bubble plot indicates the cell types that send and receive signals, and the y-axis indicates cell interactions
netVisual_bubble(cellchat, remove.isolate = FALSE, return.data = FALSE)
# Chord diagram (overall)
netVisual_chord_gene(cellchat, lab.cex = 0.8,
                     legend.pos.x = 20,
                     legend.pos.y = 20)
# Specify source and target cell groups using sources.use or targets.use
P1 <- netVisual_chord_gene(cellchat, lab.cex = 2,
                           sources.use = c("SC"),
                           legend.pos.x = 20,
                           legend.pos.y = 20)

## Bubble plot
netVisual_bubble(cellchat, targets.use = "SC", 
                 remove.isolate = TRUE) +
  theme(axis.title = element_text(size = 15),        # Axis titles
        axis.text = element_text(size = 12),         # Axis ticks
        axis.text.x = element_text(angle = 60, hjust = 0.5),  # x-axis text angle
        legend.title = element_text(size = 16, face = "bold"),  # Legend title bold
        legend.text = element_text(size = 14),                # Legend text
        legend.key.size = unit(1.5, "cm"),                    # Legend key size
        legend.key.height = unit(1, "cm"),                    # Legend key height
        legend.key.width = unit(1.2, "cm"),                    # Legend key width
        legend.spacing.y = unit(0.8, "cm"),                   # Vertical legend spacing
        legend.box.spacing = unit(0.5, "cm"))                 # Legend box spacing

# Customized plotting (if there are many ligand-receptor pairs and cell types, use remove.isolate = TRUE to remove blank rows and columns, keeping only parts with signals)
netVisual_bubble(cellchat, sources.use = "GC", 
                 targets.use = c("TC", "LC"), remove.isolate = FALSE)

# To visualize specific ligand-receptor pairs of interest, specify interaction_name from cellchat@LR$LRsig using pairLR.use (pairLR.use needs to be in dataframe format)
# Use function extractEnrichedLR to extract all important interactions (L-R pairs) and related signaling genes for a given signaling pathway (enable with geneLR.return = TRUE)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("TNF", "IL16"))
netVisual_bubble(cellchat, sources.use = c("EC", "LEC"), 
                 targets.use = c("GC", "TC", "LC", "SMC"), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE)

saveRDS(cellchat, "cell_communication.rds")




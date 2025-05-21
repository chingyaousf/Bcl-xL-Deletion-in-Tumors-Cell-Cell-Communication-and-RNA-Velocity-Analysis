# ===========================================
# Load Required Libraries
# ===========================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(CellChat)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
cellchat_path <- paste0(data_path, "cellchat_het/")

# ===========================================
# Load the Seurat Object
# ===========================================
seurat_filename <- "bclxl_annotated.rds"
bclxl <- readRDS(paste0(data_path, seurat_filename))

# Inspect the Seurat object
head(bclxl@meta.data)
print(bclxl)

# ===========================================
# Subset the Seurat Object for 'het'
# ===========================================
bclxl_het <- subset(bclxl, subset = geno == "het")

# ===========================================
# Set Up CellChat Analysis
# ===========================================
# Load CellChat database for mouse
CellChatDB <- CellChatDB.mouse

# Set cell identities and default assay
Idents(bclxl_het) <- bclxl_het$cell_type
DefaultAssay(bclxl_het) <- "RNA"

# Create CellChat object
het_cellchat <- createCellChat(object = bclxl_het, group.by = "cell_type")

# Drop unused factor levels
het_cellchat@idents <- droplevels(het_cellchat@idents)

# Set CellChat database
het_cellchat@DB <- CellChatDB

# ===========================================
# Pre-processing
# ===========================================
het_cellchat <- subsetData(het_cellchat)
het_cellchat <- identifyOverExpressedGenes(het_cellchat)
het_cellchat <- identifyOverExpressedInteractions(het_cellchat)

# Compute communication probability
het_cellchat <- computeCommunProb(het_cellchat)

# Filter communication
het_cellchat <- filterCommunication(het_cellchat, min.cells = 10)
communication.df <- subsetCommunication(het_cellchat)

# Save communication interactions
write.csv(communication.df, file = paste0(cellchat_path, "communication_interactions_het.csv"), row.names = FALSE)

# Compute pathway-level communication
het_cellchat <- computeCommunProbPathway(het_cellchat)
saveRDS(het_cellchat, file = paste0(cellchat_path, "het_cellchat.rds"))

# ===========================================
# Load Saved CellChat Object
# ===========================================
het_cellchat <- readRDS(file = paste0(cellchat_path, "het_cellchat.rds"))

# Compute aggregated network
het_cellchat <- aggregateNet(het_cellchat)

# Adjust interaction matrices
matrix <- as.data.frame(het_cellchat@net$count)
matrix[matrix == 0] <- 0.000000001  
matrix <- data.matrix(matrix)

# ===========================================
# Visualization
# ===========================================

# Interaction Number Plot
png(filename = paste0(cellchat_path, "interaction_number_het.png"), 
    width = 10, height = 12, units = "in", res = 500)
netVisual_circle(het_cellchat@net$count, top = 0.1, arrow.width = 2, arrow.size = 0.5, 
                 edge.width.max = 10, weight.scale = TRUE, label.edge = FALSE, 
                 vertex.label.cex = 1.2, title.name = "Interaction Number - het")
dev.off()

# Interaction Strength Plot
png(filename = paste0(cellchat_path, "interaction_strength_het.png"), 
    width = 10, height = 12, units = "in", res = 500)
netVisual_circle(het_cellchat@net$weight, top = 0.1, arrow.width = 2, arrow.size = 0.5, 
                 edge.width.max = 10, weight.scale = TRUE, label.edge = FALSE, 
                 vertex.label.cex = 1.2, title.name = "Interaction Strength - het")
dev.off()

# ===========================================
# Save Pathway-Level Analysis
# ===========================================
pathways.show.all <- het_cellchat@netP$pathways

for (i in 1:length(pathways.show.all)) {
  gg <- netAnalysis_contribution(het_cellchat, signaling = pathways.show.all[i])
  ggsave(paste0(cellchat_path, "LR_contribution_het/LR_contribution_het_", pathways.show.all[i], ".pdf"), 
         plot = gg, width = 13, height = 7, units = 'in', dpi = 600)
}

# ===========================================
# Violin Plots for Ligand-Receptor Expression
# ===========================================
for (i in 1:length(pathways.show.all)) {
  gg <- plotGeneExpression(het_cellchat, signaling = pathways.show.all[i], enriched.only = TRUE, type = "violin")
  ggsave(filename = paste0(cellchat_path, "violin_LR_het/violin_LR_het_", pathways.show.all[i], "_expression.png"),
         plot = gg, width = 12, height = 7, dpi = 300)
}

# ===========================================
# Network Centrality Analysis
# ===========================================
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")

# Scatter Plot for Signaling Role
png(filename = paste0(cellchat_path, "signalingRole_scatter_het.png"), 
    width = 10, height = 10, units = "in", res = 500)
netAnalysis_signalingRole_scatter(het_cellchat, title = "Signaling Role - het", label.size = 3)
dev.off()

# Heatmap for Signaling Role
png(filename = paste0(cellchat_path, "signalingRole_heatmap_het.png"), 
    width = 12, height = 12, units = "in", res = 500)
netAnalysis_signalingRole_heatmap(het_cellchat, pattern = "outgoing", font.size = 5, title = "het") +
    netAnalysis_signalingRole_heatmap(het_cellchat, pattern = "incoming", font.size = 5 , title = "het")
dev.off()

# ===========================================
# Identify Communication Patterns
# ===========================================

# Outgoing Communication Patterns
selectK(het_cellchat, pattern = "outgoing", nrun = 5)
ggsave(paste0(cellchat_path, "selectK_outgoing_het.png"))

nPatterns = 6
png(filename = paste0(cellchat_path, "outgoing_patterns_het.png"), 
    width = 10, height = 15, units = "in", res = 300)
het_cellchat <- identifyCommunicationPatterns(het_cellchat, pattern = "outgoing", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()

# Incoming Communication Patterns
selectK(het_cellchat, pattern = "incoming", nrun = 5)
ggsave(paste0(cellchat_path, "selectK_incoming_het.png"))

nPatterns = 4
png(filename = paste0(cellchat_path, "incoming_patterns_het.png"), 
    width = 10, height = 15, units = "in", res = 300)
het_cellchat <- identifyCommunicationPatterns(het_cellchat, pattern = "incoming", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()





# ===========================================
# Load Required Libraries
# ===========================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(CellChat)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
cellchat_path <- paste0(data_path, "cellchat_ko/")

# Create directory if it doesn't exist
if (!dir.exists(cellchat_path)) {
  dir.create(cellchat_path, recursive = TRUE)
}

# ===========================================
# Load the Seurat Object and Subset KO Cells
# ===========================================
seurat_filename <- "bclxl_annotated"
bclxl <- readRDS(paste0(data_path, seurat_filename, ".rds"))

# Subset KO cells
bclxl_ko <- subset(bclxl, subset = geno == "KO")

# Inspect the metadata
head(bclxl_ko@meta.data)
table(bclxl_ko$cell_type)

# ===========================================
# Set Up CellChat Object
# ===========================================
# Load CellChat database for mouse
CellChatDB <- CellChatDB.mouse

# Set cell identities and default assay
Idents(bclxl_ko) <- bclxl_ko$cell_type
DefaultAssay(bclxl_ko) <- "RNA"

# Create CellChat object
ko_cellchat <- createCellChat(object = bclxl_ko, group.by = "cell_type")

# Drop unused factor levels from cell identities
ko_cellchat@idents <- droplevels(ko_cellchat@idents)

# Check if the data matrix matches between Seurat and CellChat objects
matrix_check <- summary(ko_cellchat@data@i == bclxl_ko@assays[["RNA"]]@data@i)['FALSE']
if (is.na(matrix_check)) {
  print("Data matrices match between Seurat and CellChat.")
} else {
  print("Data matrices do not match. Check for discrepancies.")
}

# Set the CellChat database
ko_cellchat@DB <- CellChatDB

# ===========================================
# Preprocessing for CellChat
# ===========================================
ko_cellchat <- subsetData(ko_cellchat)
ko_cellchat <- identifyOverExpressedGenes(ko_cellchat)
ko_cellchat <- identifyOverExpressedInteractions(ko_cellchat)

# Compute the communication probability
ko_cellchat <- computeCommunProb(ko_cellchat)

# Infer cellular communication network
ko_cellchat <- filterCommunication(ko_cellchat, min.cells = 10)
communication.df <- subsetCommunication(ko_cellchat)

# Save the Data Frame
write.csv(communication.df, file = paste0(cellchat_path, "communication_interactions_KO.csv"), row.names = FALSE)

# Infer the cell-cell communication at a signaling pathway level
ko_cellchat <- computeCommunProbPathway(ko_cellchat)
saveRDS(ko_cellchat, file = paste0(cellchat_path, "ko_cellchat.rds"))

# ===========================================
# Load the saved CellChat object
# ===========================================
ko_cellchat <- readRDS(file = paste0(cellchat_path, "ko_cellchat.rds"))

# ===========================================
# Compute and Visualize Cell-Cell Communication Network
# ===========================================
ko_cellchat <- aggregateNet(ko_cellchat)

# Adjust the interaction matrices and proceed with visualizations
matrix <- as.data.frame(ko_cellchat@net$count)
matrix[matrix == 0] <- 0.000000001  # Replace zeros with a very small value
matrix <- data.matrix(matrix)

# Save adjusted matrix if needed
# write.csv(matrix, file = paste0(cellchat_path, "interaction_matrix_KO.csv"), row.names = FALSE)

# Visualize the number of interactions with weights/strengths
png(filename = paste0(cellchat_path, "interaction_number_KO.png"), 
    width = 10, height = 12, units = "in", res = 500)
netVisual_circle(ko_cellchat@net$count, top = 0.1, arrow.width = 2, arrow.size = 0.5, 
                 edge.width.max = 10, weight.scale = TRUE, label.edge = FALSE, 
                 vertex.label.cex = 1.2, title.name = "Interaction Number - KO")
dev.off()

png(filename = paste0(cellchat_path, "interaction_strength_KO.png"), 
    width = 10, height = 12, units = "in", res = 500)
netVisual_circle(ko_cellchat@net$weight, top = 0.1, arrow.width = 2, arrow.size = 0.5, 
                 edge.width.max = 10, weight.scale = TRUE, label.edge = FALSE, 
                 vertex.label.cex = 1.2, title.name = "Interaction Strength - KO")
dev.off()

# ===========================================
# Save Plots for Each Inferred Network
# ===========================================
pathways.show.all <- ko_cellchat@netP$pathways
for (i in 1:length(pathways.show.all)) {
  gg <- netAnalysis_contribution(ko_cellchat, signaling = pathways.show.all[i])
  ggsave(paste0(cellchat_path, "LR_contribution_KO/LR_contribution_KO_", pathways.show.all[i], ".pdf"), 
         plot = gg, width = 13, height = 7, units = 'in', dpi = 600)
}

# ===========================================
# Generate Violin Plots for Ligand-Receptor Gene Expression
# ===========================================
for (i in 1:length(pathways.show.all)) {
  # Generate violin plot
  gg <- plotGeneExpression(ko_cellchat, signaling = pathways.show.all[i], enriched.only = TRUE, type = "violin")

  # Save violin plot
  ggsave(filename = paste0(cellchat_path, "violin_LR_KO/violin_LR_KO_", pathways.show.all[i], ".png"),
         plot = gg, width = 12, height = 7, dpi = 300)
}

# ===========================================
# Compute and Visualize Network Centrality Scores
# ===========================================
ko_cellchat <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# Scatter plot for signaling role analysis
png(filename = paste0(cellchat_path, "signalingRole_scatter_KO.png"), 
    width = 5, height = 5, units = "in", res = 500)
netAnalysis_signalingRole_scatter(ko_cellchat, title = "Signaling Role - KO")
dev.off()

# Heatmap for signaling role analysis
png(filename = paste0(cellchat_path, "signalingRole_heatmap_KO.png"), 
    width = 12, height = 12, units = "in", res = 500)
netAnalysis_signalingRole_heatmap(ko_cellchat, pattern = "outgoing", font.size = 5, title = "KO") +
netAnalysis_signalingRole_heatmap(ko_cellchat, pattern = "incoming", font.size = 5 , title = "KO")
dev.off()

# ===========================================
# Identify Outgoing and Incoming Communication Patterns
# ===========================================
selectK(ko_cellchat, pattern = "outgoing", nrun = 5)
ggsave(paste0(cellchat_path, "selectK_outgoing_KO.png"))
nPatterns = 6
png(filename = paste0(cellchat_path, "outgoing_patterns_KO.png"), 
    width = 10, height = 15, units = "in", res = 300)
ko_cellchat <- identifyCommunicationPatterns(ko_cellchat, pattern = "outgoing", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()

selectK(ko_cellchat, pattern = "incoming", nrun = 5)
ggsave(paste0(cellchat_path, "selectK_incoming_KO.png"))
nPatterns = 4
png(filename = paste0(cellchat_path, "incoming_patterns_KO.png"), 
    width = 10, height = 15, units = "in", res = 300)
ko_cellchat <- identifyCommunicationPatterns(ko_cellchat, pattern = "incoming", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()




# ===========================================
# Load Required Libraries
# ===========================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(CellChat)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
cellchat_path <- paste0(data_path, "cellchat_het/")

# ===========================================
# Load the Seurat Object
# ===========================================
seurat_filename <- "bclxl_annotated.rds"
bclxl <- readRDS(paste0(data_path, seurat_filename))

# Inspect the Seurat object
head(bclxl@meta.data)
print(bclxl)

# ===========================================
# Subset the Seurat Object for 'het'
# ===========================================
bclxl_het <- subset(bclxl, subset = geno == "het")

# ===========================================
# Set Up CellChat Analysis
# ===========================================
# Load CellChat database for mouse
CellChatDB <- CellChatDB.mouse

# Set cell identities and default assay
Idents(bclxl_het) <- bclxl_het$cell_type
DefaultAssay(bclxl_het) <- "RNA"

# Create CellChat object
het_cellchat <- createCellChat(object = bclxl_het, group.by = "cell_type")

# Drop unused factor levels
het_cellchat@idents <- droplevels(het_cellchat@idents)

# Set CellChat database
het_cellchat@DB <- CellChatDB

# ===========================================
# Pre-processing
# ===========================================
het_cellchat <- subsetData(het_cellchat)
het_cellchat <- identifyOverExpressedGenes(het_cellchat)
het_cellchat <- identifyOverExpressedInteractions(het_cellchat)

# Compute communication probability
het_cellchat <- computeCommunProb(het_cellchat)

# Filter communication
het_cellchat <- filterCommunication(het_cellchat, min.cells = 10)
communication.df <- subsetCommunication(het_cellchat)

# Save communication interactions
write.csv(communication.df, file = paste0(cellchat_path, "communication_interactions_het.csv"), row.names = FALSE)

# Compute pathway-level communication
het_cellchat <- computeCommunProbPathway(het_cellchat)
saveRDS(het_cellchat, file = paste0(cellchat_path, "het_cellchat.rds"))

# ===========================================
# Load Saved CellChat Object
# ===========================================
het_cellchat <- readRDS(file = paste0(cellchat_path, "het_cellchat.rds"))

# Compute aggregated network
het_cellchat <- aggregateNet(het_cellchat)

# Adjust interaction matrices
matrix <- as.data.frame(het_cellchat@net$count)
matrix[matrix == 0] <- 0.000000001  
matrix <- data.matrix(matrix)

# ===========================================
# Visualization
# ===========================================

# Interaction Number Plot
png(filename = paste0(cellchat_path, "interaction_number_het.png"), 
    width = 10, height = 12, units = "in", res = 500)
netVisual_circle(het_cellchat@net$count, top = 0.1, arrow.width = 2, arrow.size = 0.5, 
                 edge.width.max = 10, weight.scale = TRUE, label.edge = FALSE, 
                 vertex.label.cex = 1.2, title.name = "Interaction Number - het")
dev.off()

# Interaction Strength Plot
png(filename = paste0(cellchat_path, "interaction_strength_het.png"), 
    width = 10, height = 12, units = "in", res = 500)
netVisual_circle(het_cellchat@net$weight, top = 0.1, arrow.width = 2, arrow.size = 0.5, 
                 edge.width.max = 10, weight.scale = TRUE, label.edge = FALSE, 
                 vertex.label.cex = 1.2, title.name = "Interaction Strength - het")
dev.off()

# ===========================================
# Save Pathway-Level Analysis
# ===========================================
pathways.show.all <- het_cellchat@netP$pathways

for (i in 1:length(pathways.show.all)) {
  gg <- netAnalysis_contribution(het_cellchat, signaling = pathways.show.all[i])
  ggsave(paste0(cellchat_path, "LR_contribution_het/LR_contribution_het_", pathways.show.all[i], ".pdf"), 
         plot = gg, width = 13, height = 7, units = 'in', dpi = 600)
}

# ===========================================
# Violin Plots for Ligand-Receptor Expression
# ===========================================
for (i in 1:length(pathways.show.all)) {
  gg <- plotGeneExpression(het_cellchat, signaling = pathways.show.all[i], enriched.only = TRUE, type = "violin")
  ggsave(filename = paste0(cellchat_path, "violin_LR_het/violin_LR_het_", pathways.show.all[i], "_expression.png"),
         plot = gg, width = 12, height = 7, dpi = 300)
}

# ===========================================
# Network Centrality Analysis
# ===========================================
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")

# Scatter Plot for Signaling Role
png(filename = paste0(cellchat_path, "signalingRole_scatter_het.png"), 
    width = 10, height = 10, units = "in", res = 500)
netAnalysis_signalingRole_scatter(het_cellchat, title = "Signaling Role - het", label.size = 3)
dev.off()

# Heatmap for Signaling Role
png(filename = paste0(cellchat_path, "signalingRole_heatmap_het.png"), 
    width = 12, height = 12, units = "in", res = 500)
netAnalysis_signalingRole_heatmap(het_cellchat, pattern = "outgoing", font.size = 5, title = "het") +
    netAnalysis_signalingRole_heatmap(het_cellchat, pattern = "incoming", font.size = 5 , title = "het")
dev.off()

# ===========================================
# Identify Communication Patterns
# ===========================================

# Outgoing Communication Patterns
selectK(het_cellchat, pattern = "outgoing", nrun = 5)
ggsave(paste0(cellchat_path, "selectK_outgoing_het.png"))

nPatterns = 6
png(filename = paste0(cellchat_path, "outgoing_patterns_het.png"), 
    width = 10, height = 15, units = "in", res = 300)
het_cellchat <- identifyCommunicationPatterns(het_cellchat, pattern = "outgoing", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()

# Incoming Communication Patterns
selectK(het_cellchat, pattern = "incoming", nrun = 5)
ggsave(paste0(cellchat_path, "selectK_incoming_het.png"))

nPatterns = 4
png(filename = paste0(cellchat_path, "incoming_patterns_het.png"), 
    width = 10, height = 15, units = "in", res = 300)
het_cellchat <- identifyCommunicationPatterns(het_cellchat, pattern = "incoming", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()




### bclxl_annotated_noRBC, remove Red_Blood_Cells, cellchat_ko_02
# ===========================================
# Load Required Libraries
# ===========================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(CellChat)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
cellchat_path <- paste0(data_path, "cellchat_ko_02/")

# Create directory if it doesn't exist
if (!dir.exists(cellchat_path)) {
  dir.create(cellchat_path, recursive = TRUE)
}

# ===========================================
# Load the Seurat Object and Subset KO Cells
# ===========================================
seurat_filename <- "bclxl_annotated_noRBC"
bclxl <- readRDS(paste0(data_path, seurat_filename, ".rds"))

# Subset KO cells
bclxl_ko <- subset(bclxl, subset = geno == "KO")

# Inspect the metadata
head(bclxl_ko@meta.data)
table(bclxl_ko$cell_type)

# ===========================================
# Set Up CellChat Object
# ===========================================
# Load CellChat database for mouse
CellChatDB <- CellChatDB.mouse

# Set cell identities and default assay
Idents(bclxl_ko) <- bclxl_ko$cell_type
DefaultAssay(bclxl_ko) <- "RNA"

# Create CellChat object
ko_cellchat <- createCellChat(object = bclxl_ko, group.by = "cell_type")

# Drop unused factor levels from cell identities
ko_cellchat@idents <- droplevels(ko_cellchat@idents)

# Check if the data matrix matches between Seurat and CellChat objects
matrix_check <- summary(ko_cellchat@data@i == bclxl_ko@assays[["RNA"]]@data@i)['FALSE']
if (is.na(matrix_check)) {
  print("Data matrices match between Seurat and CellChat.")
} else {
  print("Data matrices do not match. Check for discrepancies.")
}

# Set the CellChat database
ko_cellchat@DB <- CellChatDB

# ===========================================
# Preprocessing for CellChat
# ===========================================
ko_cellchat <- subsetData(ko_cellchat)
ko_cellchat <- identifyOverExpressedGenes(ko_cellchat)
ko_cellchat <- identifyOverExpressedInteractions(ko_cellchat)

# Compute the communication probability
ko_cellchat <- computeCommunProb(ko_cellchat)

# Infer cellular communication network
ko_cellchat <- filterCommunication(ko_cellchat, min.cells = 10)
communication.df <- subsetCommunication(ko_cellchat)

# Save the Data Frame
write.csv(communication.df, file = paste0(cellchat_path, "communication_interactions_KO.csv"), row.names = FALSE)

# Infer the cell-cell communication at a signaling pathway level
ko_cellchat <- computeCommunProbPathway(ko_cellchat)
saveRDS(ko_cellchat, file = paste0(cellchat_path, "ko_cellchat_02.rds"))

# ===========================================
# Load the saved CellChat object
# ===========================================
ko_cellchat <- readRDS(file = paste0(cellchat_path, "ko_cellchat_02.rds"))

# ===========================================
# Compute and Visualize Cell-Cell Communication Network
# ===========================================
ko_cellchat <- aggregateNet(ko_cellchat)

# Adjust the interaction matrices and proceed with visualizations
matrix <- as.data.frame(ko_cellchat@net$count)
matrix[matrix == 0] <- 0.000000001  # Replace zeros with a very small value
matrix <- data.matrix(matrix)

# Save adjusted matrix if needed
# write.csv(matrix, file = paste0(cellchat_path, "interaction_matrix_KO.csv"), row.names = FALSE)

# Visualize the number of interactions with weights/strengths
png(filename = paste0(cellchat_path, "interaction_number_KO_02.png"), 
    width = 13, height = 12, units = "in", res = 500)
netVisual_circle(ko_cellchat@net$count, top = 0.1, arrow.width = 2, arrow.size = 0.5, 
                 edge.width.max = 10, weight.scale = TRUE, label.edge = FALSE, 
                 vertex.label.cex = 1.9, title.name = "Interaction Number - KO")
dev.off()

png(filename = paste0(cellchat_path, "interaction_strength_KO_02.png"), 
    width = 13, height = 12, units = "in", res = 500)
netVisual_circle(ko_cellchat@net$weight, top = 0.1, arrow.width = 2, arrow.size = 0.5, 
                 edge.width.max = 10, weight.scale = TRUE, label.edge = FALSE, 
                 vertex.label.cex = 1.9, title.name = "Interaction Strength - KO")
dev.off()

# ===========================================
# Save Plots for Each Inferred Network
# ===========================================
pathways.show.all <- ko_cellchat@netP$pathways
for (i in 1:length(pathways.show.all)) {
  gg <- netAnalysis_contribution(ko_cellchat, signaling = pathways.show.all[i])
  ggsave(paste0(cellchat_path, "LR_contribution_KO/LR_contribution_KO_", pathways.show.all[i], ".pdf"), 
         plot = gg, width = 13, height = 7, units = 'in', dpi = 600)
}

# ===========================================
# Generate Violin Plots for Ligand-Receptor Gene Expression
# ===========================================
for (i in 1:length(pathways.show.all)) {
  # Generate violin plot
  gg <- plotGeneExpression(ko_cellchat, signaling = pathways.show.all[i], enriched.only = TRUE, type = "violin")

  # Save violin plot
  ggsave(filename = paste0(cellchat_path, "violin_LR_KO/violin_LR_KO_", pathways.show.all[i], ".png"),
         plot = gg, width = 12, height = 7, dpi = 300)
}

# ===========================================
# Compute and Visualize Network Centrality Scores
# ===========================================
ko_cellchat <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# Scatter plot for signaling role analysis
png(filename = paste0(cellchat_path, "signalingRole_scatter_KO.png"), 
    width = 5, height = 5, units = "in", res = 500)
netAnalysis_signalingRole_scatter(ko_cellchat, title = "Signaling Role - KO")
dev.off()

# Heatmap for signaling role analysis
png(filename = paste0(cellchat_path, "signalingRole_heatmap_KO.png"), 
    width = 12, height = 12, units = "in", res = 500)
netAnalysis_signalingRole_heatmap(ko_cellchat, pattern = "outgoing", font.size = 5, title = "KO") +
netAnalysis_signalingRole_heatmap(ko_cellchat, pattern = "incoming", font.size = 5 , title = "KO")
dev.off()

# ===========================================
# Identify Outgoing and Incoming Communication Patterns
# ===========================================
selectK(ko_cellchat, pattern = "outgoing", nrun = 5)
ggsave(paste0(cellchat_path, "selectK_outgoing_KO.png"))
nPatterns = 6
png(filename = paste0(cellchat_path, "outgoing_patterns_KO.png"), 
    width = 10, height = 15, units = "in", res = 300)
ko_cellchat <- identifyCommunicationPatterns(ko_cellchat, pattern = "outgoing", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()

selectK(ko_cellchat, pattern = "incoming", nrun = 5)
ggsave(paste0(cellchat_path, "selectK_incoming_KO.png"))
nPatterns = 4
png(filename = paste0(cellchat_path, "incoming_patterns_KO.png"), 
    width = 10, height = 15, units = "in", res = 300)
ko_cellchat <- identifyCommunicationPatterns(ko_cellchat, pattern = "incoming", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()




### bclxl_annotated_noRBC, remove Red_Blood_Cells, cellchat_het_02
# ===========================================
# Load Required Libraries
# ===========================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(CellChat)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
cellchat_path <- paste0(data_path, "cellchat_het_02/")

# ===========================================
# Load the Seurat Object
# ===========================================
seurat_filename <- "bclxl_annotated_noRBC.rds"
bclxl <- readRDS(paste0(data_path, seurat_filename))

# Inspect the Seurat object
head(bclxl@meta.data)
print(bclxl)

# ===========================================
# Subset the Seurat Object for 'het'
# ===========================================
bclxl_het <- subset(bclxl, subset = geno == "het")

# ===========================================
# Set Up CellChat Analysis
# ===========================================
# Load CellChat database for mouse
CellChatDB <- CellChatDB.mouse

# Set cell identities and default assay
Idents(bclxl_het) <- bclxl_het$cell_type
DefaultAssay(bclxl_het) <- "RNA"

# Create CellChat object
het_cellchat <- createCellChat(object = bclxl_het, group.by = "cell_type")

# Drop unused factor levels
het_cellchat@idents <- droplevels(het_cellchat@idents)

# Set CellChat database
het_cellchat@DB <- CellChatDB

# ===========================================
# Pre-processing
# ===========================================
het_cellchat <- subsetData(het_cellchat)
het_cellchat <- identifyOverExpressedGenes(het_cellchat)
het_cellchat <- identifyOverExpressedInteractions(het_cellchat)

# Compute communication probability
het_cellchat <- computeCommunProb(het_cellchat)

# Filter communication
het_cellchat <- filterCommunication(het_cellchat, min.cells = 10)
communication.df <- subsetCommunication(het_cellchat)

# Save communication interactions
write.csv(communication.df, file = paste0(cellchat_path, "communication_interactions_het.csv"), row.names = FALSE)

# Compute pathway-level communication
het_cellchat <- computeCommunProbPathway(het_cellchat)
saveRDS(het_cellchat, file = paste0(cellchat_path, "het_cellchat_02.rds"))

# ===========================================
# Load Saved CellChat Object
# ===========================================
het_cellchat <- readRDS(file = paste0(cellchat_path, "het_cellchat_02.rds"))

# Compute aggregated network
het_cellchat <- aggregateNet(het_cellchat)

# Adjust interaction matrices
matrix <- as.data.frame(het_cellchat@net$count)
matrix[matrix == 0] <- 0.000000001  
matrix <- data.matrix(matrix)

# ===========================================
# Visualization
# ===========================================

# Interaction Number Plot
png(filename = paste0(cellchat_path, "interaction_number_het_02.png"), 
    width = 13, height = 12, units = "in", res = 500)
netVisual_circle(het_cellchat@net$count, top = 0.1, arrow.width = 2, arrow.size = 0.5, 
                 edge.width.max = 10, weight.scale = TRUE, label.edge = FALSE, 
                 vertex.label.cex = 1.9, title.name = "Interaction Number - het")
dev.off()

# Interaction Strength Plot
png(filename = paste0(cellchat_path, "interaction_strength_het_02.png"), 
    width = 13, height = 12, units = "in", res = 500)
netVisual_circle(het_cellchat@net$weight, top = 0.1, arrow.width = 2, arrow.size = 0.5, 
                 edge.width.max = 10, weight.scale = TRUE, label.edge = FALSE, 
                 vertex.label.cex = 1.9, title.name = "Interaction Strength - het")
dev.off()

# ===========================================
# Save Pathway-Level Analysis
# ===========================================
pathways.show.all <- het_cellchat@netP$pathways

for (i in 1:length(pathways.show.all)) {
  gg <- netAnalysis_contribution(het_cellchat, signaling = pathways.show.all[i])
  ggsave(paste0(cellchat_path, "LR_contribution_het/LR_contribution_het_", pathways.show.all[i], ".pdf"), 
         plot = gg, width = 13, height = 7, units = 'in', dpi = 600)
}

# ===========================================
# Violin Plots for Ligand-Receptor Expression
# ===========================================
for (i in 1:length(pathways.show.all)) {
  gg <- plotGeneExpression(het_cellchat, signaling = pathways.show.all[i], enriched.only = TRUE, type = "violin")
  ggsave(filename = paste0(cellchat_path, "violin_LR_het/violin_LR_het_", pathways.show.all[i], "_expression.png"),
         plot = gg, width = 12, height = 7, dpi = 300)
}

# ===========================================
# Network Centrality Analysis
# ===========================================
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")

# Scatter Plot for Signaling Role
png(filename = paste0(cellchat_path, "signalingRole_scatter_het.png"), 
    width = 10, height = 10, units = "in", res = 500)
netAnalysis_signalingRole_scatter(het_cellchat, title = "Signaling Role - het", label.size = 3)
dev.off()

# Heatmap for Signaling Role
png(filename = paste0(cellchat_path, "signalingRole_heatmap_het.png"), 
    width = 12, height = 12, units = "in", res = 500)
netAnalysis_signalingRole_heatmap(het_cellchat, pattern = "outgoing", font.size = 5, title = "het") +
    netAnalysis_signalingRole_heatmap(het_cellchat, pattern = "incoming", font.size = 5 , title = "het")
dev.off()

# ===========================================
# Identify Communication Patterns
# ===========================================

# Outgoing Communication Patterns
selectK(het_cellchat, pattern = "outgoing", nrun = 5)
ggsave(paste0(cellchat_path, "selectK_outgoing_het.png"))

nPatterns = 6
png(filename = paste0(cellchat_path, "outgoing_patterns_het.png"), 
    width = 10, height = 15, units = "in", res = 300)
het_cellchat <- identifyCommunicationPatterns(het_cellchat, pattern = "outgoing", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()

# Incoming Communication Patterns
selectK(het_cellchat, pattern = "incoming", nrun = 5)
ggsave(paste0(cellchat_path, "selectK_incoming_het.png"))

nPatterns = 4
png(filename = paste0(cellchat_path, "incoming_patterns_het.png"), 
    width = 10, height = 15, units = "in", res = 300)
het_cellchat <- identifyCommunicationPatterns(het_cellchat, pattern = "incoming", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()

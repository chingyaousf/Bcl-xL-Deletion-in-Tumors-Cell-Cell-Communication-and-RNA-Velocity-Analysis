
# ==============================================================================
# APP_circle_comparison between Het and KO for All Cell Type Interactions
# ==============================================================================
# Title: Generate Chord Diagrams for the APP Signaling Pathway
# Description: This script generates chord diagrams for the APP signaling pathway 
# to visualize the interactions between all cell types in both Het and KO datasets.
# No specific cell type is highlighted.
# ==============================================================================

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko/")


# ===========================================
# Load the Saved CellChat Objects
# ===========================================
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het/het_cellchat.rds"))
ko_cellchat <- readRDS(file = paste0(data_path, "cellchat_ko/ko_cellchat.rds"))

# ===========================================
# Aggregate Communication Networks for Each Dataset
# ===========================================
het_cellchat <- aggregateNet(het_cellchat)
ko_cellchat <- aggregateNet(ko_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Het = het_cellchat, KO = ko_cellchat)

# ===========================================
# Define the Signaling Pathway to Visualize
# ===========================================
pathways.show <- c("APP")  # APP pathway for visualization

# No specific cell type is highlighted
group.cellType <- rep("All", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)

# ---- STEP 1: Display Chord Diagrams Interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Het and KO
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

# ---- STEP 2: Save Chord Diagrams after Review ----
# Save the Chord Diagrams to a PNG file
png(filename = paste0(comparison_path, "APP_chord_comparison_het_ko.png"),
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Het and KO
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

dev.off()





# ==============================================================================
# Generate Chord Diagrams for CADM, PCDH, and PSAP Signaling Pathways
# Comparing Het and KO in CellChat
# ==============================================================================

# Load required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko/")

# Create directory if it doesn't exist
if (!dir.exists(comparison_path)) {
  dir.create(comparison_path, recursive = TRUE)
}

# ===========================================
# Load the Saved CellChat Objects
# ===========================================
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het/het_cellchat.rds"))
ko_cellchat <- readRDS(file = paste0(data_path, "cellchat_ko/ko_cellchat.rds"))

# ===========================================
# Aggregate Communication Networks for Each Dataset
# ===========================================
het_cellchat <- aggregateNet(het_cellchat)
ko_cellchat <- aggregateNet(ko_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Het = het_cellchat, KO = ko_cellchat)

# ===========================================
# Define the Signaling Pathways to Visualize
# ===========================================
pathways.show <- c("CADM", "PCDH", "PSAP")  # List of signaling pathways

# No specific cell type is highlighted (treat all equally)
group.cellType <- rep("All", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)

# ---- STEP 1: Display Chord Diagrams Interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 3), xpd = TRUE)  # 3 plots for CADM, PCDH, and PSAP
for (p in pathways.show) {
  for (i in 1:length(object.list)) {
    netVisual_chord_cell(object.list[[i]], 
                         signaling = p, 
                         group = group.cellType, 
                         title.name = paste0(p, " signaling network - ", names(object.list)[i]))
  }
}

# ---- STEP 2: Save Chord Diagrams after Review ----
# Save each Chord Diagram separately as PNG files
for (p in pathways.show) {
  png(filename = paste0(comparison_path, p, "_chord_comparison_het_ko.png"),
      width = 16, height = 8, units = "in", res = 300)

  # Re-create the Chord Diagrams for saving
  par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Het and KO
  for (i in 1:length(object.list)) {
    netVisual_chord_cell(object.list[[i]], 
                         signaling = p, 
                         group = group.cellType, 
                         title.name = paste0(p, " signaling network - ", names(object.list)[i]))
  }

  dev.off()
}





# ==============================================================================
# Generate Chord Diagrams for MPZ, Cholesterol, SPP1, CLDN, and VTN
# Comparing Het and KO in CellChat
# ==============================================================================

# Load required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko/")


# ===========================================
# Load the Saved CellChat Objects
# ===========================================
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het/het_cellchat.rds"))
ko_cellchat <- readRDS(file = paste0(data_path, "cellchat_ko/ko_cellchat.rds"))

# ===========================================
# Aggregate Communication Networks for Each Dataset
# ===========================================
het_cellchat <- aggregateNet(het_cellchat)
ko_cellchat <- aggregateNet(ko_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Het = het_cellchat, KO = ko_cellchat)

# ===========================================
# Define the Signaling Pathways to Visualize
# ===========================================
pathways.show <- c("MPZ", "Cholesterol", "SPP1", "CLDN", "VTN")  # Updated pathways

# No specific cell type is highlighted (treat all equally)
group.cellType <- rep("All", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)

# ---- STEP 1: Display Chord Diagrams Interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 3), xpd = TRUE)  # Show 3 plots at a time for clarity

for (p in pathways.show) {
  for (i in 1:length(object.list)) {
    if (p %in% object.list[[i]]@netP$pathways) {  # Skip if pathway is not significant
      netVisual_chord_cell(object.list[[i]], 
                           signaling = p, 
                           group = group.cellType, 
                           title.name = paste0(p, " signaling network - ", names(object.list)[i]))
    } else {
      message(paste("Skipping", p, "for", names(object.list)[i], "due to no significant communication."))
    }
  }
}

# ---- STEP 2: Save Chord Diagrams after Review ----
# Save each Chord Diagram separately as PNG files
for (p in pathways.show) {
  if (any(p %in% het_cellchat@netP$pathways, p %in% ko_cellchat@netP$pathways)) {  # Check if any dataset has significant communication
    png(filename = paste0(comparison_path, p, "_chord_comparison_het_ko.png"),
        width = 16, height = 8, units = "in", res = 300)

    # Re-create the Chord Diagrams for saving
    par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Het and KO
    for (i in 1:length(object.list)) {
      if (p %in% object.list[[i]]@netP$pathways) {
        netVisual_chord_cell(object.list[[i]], 
                             signaling = p, 
                             group = group.cellType, 
                             title.name = paste0(p, " signaling network - ", names(object.list)[i]))
      } else {
        message(paste("Skipping", p, "for", names(object.list)[i], "due to no significant communication."))
      }
    }

    dev.off()
  } else {
    message(paste("Skipping", p, "entirely due to no significant communication in both Het and KO."))
  }
}




#### remove Red_Blood_Cells
# ==============================================================================
# APP_circle_comparison between Het and KO for All Cell Type Interactions
# ==============================================================================
# Title: Generate Chord Diagrams for the APP Signaling Pathway
# Description: This script generates chord diagrams for the APP signaling pathway 
# to visualize the interactions between all cell types in both Het and KO datasets.
# No specific cell type is highlighted.
# ==============================================================================

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko_02/")


# ===========================================
# Load the Saved CellChat Objects
# ===========================================
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het_02/het_cellchat_02.rds"))
ko_cellchat <- readRDS(file = paste0(data_path, "cellchat_ko_02/ko_cellchat_02.rds"))

# ===========================================
# Aggregate Communication Networks for Each Dataset
# ===========================================
het_cellchat <- aggregateNet(het_cellchat)
ko_cellchat <- aggregateNet(ko_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Het = het_cellchat, KO = ko_cellchat)

# ===========================================
# Define the Signaling Pathway to Visualize
# ===========================================
pathways.show <- c("APP")  # APP pathway for visualization

# No specific cell type is highlighted
group.cellType <- rep("All", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)

# ---- STEP 1: Display Chord Diagrams Interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Het and KO
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

# ---- STEP 2: Save Chord Diagrams after Review ----
# Save the Chord Diagrams to a PNG file
png(filename = paste0(comparison_path, "APP_chord_comparison_het_ko.png"),
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Het and KO
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

dev.off()





# ==============================================================================
# Generate Chord Diagrams for CADM, PCDH, and PSAP Signaling Pathways
# Comparing Het and KO in CellChat
# ==============================================================================

# Load required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko_02/")

# Create directory if it doesn't exist
if (!dir.exists(comparison_path)) {
  dir.create(comparison_path, recursive = TRUE)
}

# ===========================================
# Load the Saved CellChat Objects
# ===========================================
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het_02/het_cellchat_02.rds"))
ko_cellchat <- readRDS(file = paste0(data_path, "cellchat_ko_02/ko_cellchat_02.rds"))

# ===========================================
# Aggregate Communication Networks for Each Dataset
# ===========================================
het_cellchat <- aggregateNet(het_cellchat)
ko_cellchat <- aggregateNet(ko_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Het = het_cellchat, KO = ko_cellchat)

# ===========================================
# Define the Signaling Pathways to Visualize
# ===========================================
pathways.show <- c("PCDH","ESAM","GAP","CNTN")  # List of signaling pathways

# No specific cell type is highlighted (treat all equally)
group.cellType <- rep("All", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)

# ---- STEP 1: Display Chord Diagrams Interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 3), xpd = TRUE)  # 3 plots for CADM, PCDH, and PSAP
for (p in pathways.show) {
  for (i in 1:length(object.list)) {
    netVisual_chord_cell(object.list[[i]], 
                         signaling = p, 
                         group = group.cellType, 
                         title.name = paste0(p, " signaling network - ", names(object.list)[i]))
  }
}

# ---- STEP 2: Save Chord Diagrams after Review ----
# Save each Chord Diagram separately as PNG files
for (p in pathways.show) {
  png(filename = paste0(comparison_path, p, "_chord_comparison_het_ko.png"),
      width = 16, height = 8, units = "in", res = 300)

  # Re-create the Chord Diagrams for saving
  par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Het and KO
  for (i in 1:length(object.list)) {
    netVisual_chord_cell(object.list[[i]], 
                         signaling = p, 
                         group = group.cellType, 
                         title.name = paste0(p, " signaling network - ", names(object.list)[i]))
  }

  dev.off()
}





# ==============================================================================
# Generate Chord Diagrams for MPZ, Cholesterol, SPP1, CLDN, and VTN
# Comparing Het and KO in CellChat
# ==============================================================================

# Load required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko_02/")


# ===========================================
# Load the Saved CellChat Objects
# ===========================================
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het_02/het_cellchat_02.rds"))
ko_cellchat <- readRDS(file = paste0(data_path, "cellchat_ko_02/ko_cellchat_02.rds"))

# ===========================================
# Aggregate Communication Networks for Each Dataset
# ===========================================
het_cellchat <- aggregateNet(het_cellchat)
ko_cellchat <- aggregateNet(ko_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Het = het_cellchat, KO = ko_cellchat)

# ===========================================
# Define the Signaling Pathways to Visualize
# ===========================================
pathways.show <- c("MPZ", "Cholesterol", "SPP1", "CLDN", "VTN")  # Updated pathways

# No specific cell type is highlighted (treat all equally)
group.cellType <- rep("All", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)

# ---- STEP 1: Display Chord Diagrams Interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 3), xpd = TRUE)  # Show 3 plots at a time for clarity

for (p in pathways.show) {
  for (i in 1:length(object.list)) {
    if (p %in% object.list[[i]]@netP$pathways) {  # Skip if pathway is not significant
      netVisual_chord_cell(object.list[[i]], 
                           signaling = p, 
                           group = group.cellType, 
                           title.name = paste0(p, " signaling network - ", names(object.list)[i]))
    } else {
      message(paste("Skipping", p, "for", names(object.list)[i], "due to no significant communication."))
    }
  }
}

# ---- STEP 2: Save Chord Diagrams after Review ----
# Save each Chord Diagram separately as PNG files
for (p in pathways.show) {
  if (any(p %in% het_cellchat@netP$pathways, p %in% ko_cellchat@netP$pathways)) {  # Check if any dataset has significant communication
    png(filename = paste0(comparison_path, p, "_chord_comparison_het_ko.png"),
        width = 16, height = 8, units = "in", res = 300)

    # Re-create the Chord Diagrams for saving
    par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Het and KO
    for (i in 1:length(object.list)) {
      if (p %in% object.list[[i]]@netP$pathways) {
        netVisual_chord_cell(object.list[[i]], 
                             signaling = p, 
                             group = group.cellType, 
                             title.name = paste0(p, " signaling network - ", names(object.list)[i]))
      } else {
        message(paste("Skipping", p, "for", names(object.list)[i], "due to no significant communication."))
      }
    }

    dev.off()
  } else {
    message(paste("Skipping", p, "entirely due to no significant communication in both Het and KO."))
  }
}



# ==============================================================================
# Generate Chord Diagrams for "PCDH","ESAM","CNTN","CDH5","CADM","NCAM","GAP","FLRT","ApoE","APP","PSAP","PTN","MK","PECAM1","EDN","CypA" Signaling Pathways
# Comparing Het and KO in CellChat
# ==============================================================================

# Load required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko_02/")

# Create directory if it doesn't exist
if (!dir.exists(comparison_path)) {
  dir.create(comparison_path, recursive = TRUE)
}

# ===========================================
# Load the Saved CellChat Objects
# ===========================================
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het_02/het_cellchat_02.rds"))
ko_cellchat <- readRDS(file = paste0(data_path, "cellchat_ko_02/ko_cellchat_02.rds"))

# ===========================================
# Aggregate Communication Networks for Each Dataset
# ===========================================
het_cellchat <- aggregateNet(het_cellchat)
ko_cellchat <- aggregateNet(ko_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Het = het_cellchat, KO = ko_cellchat)

# ===========================================
# Define the Signaling Pathways to Visualize
# ===========================================
#pathways.show <- c("PCDH","ESAM","CNTN","CDH5","CADM","NCAM","GAP","FLRT","ApoE","APP","PSAP","PTN","MK","PECAM1","EDN","CypA")  # List of signaling pathways
pathways.show <- c("PCDH")  # List of signaling pathways

# No specific cell type is highlighted (treat all equally)
group.cellType <- rep("All", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)

# ---- STEP 1: Display Chord Diagrams Interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 3), xpd = TRUE)  # 3 plots for CADM, PCDH, and PSAP
for (p in pathways.show) {
  for (i in 1:length(object.list)) {
    netVisual_chord_cell(object.list[[i]], 
                         signaling = p, 
                         group = group.cellType, 
                         title.name = paste0(p, " signaling network - ", names(object.list)[i]))
  }
}

# ---- STEP 2: Save Chord Diagrams after Review ----
# Save each Chord Diagram separately as PNG files
for (p in pathways.show) {
  png(filename = paste0(comparison_path, p, "_chord_comparison_het_ko_02.png"),
      width = 16, height = 8, units = "in", res = 300)



  # Re-create the Chord Diagrams for saving
   par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Het and KO
      
  for (i in 1:length(object.list)) {
    netVisual_chord_cell(object.list[[i]],
                         signaling = p, 
                         group = group.cellType, 
                         title.name = paste0(p, " signaling network - ", names(object.list)[i]))
  }

  dev.off()
}

# ===========================================
# Load Required Libraries
# ===========================================
rm(list = ls())

library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(dplyr)

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

# Check Data Consistency: Ensure both have sufficient overlap in signaling pathways and cell types
signaling_het <- het_cellchat@netP$pathways
signaling_ko <- ko_cellchat@netP$pathways

# Prepare the CellChat objects by aggregating the communication networks
het_cellchat <- aggregateNet(het_cellchat)
ko_cellchat <- aggregateNet(ko_cellchat)

# Calculate the centrality scores for the signaling roles
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# Create a list of CellChat objects for comparison
object.list <- list(Het = het_cellchat, KO = ko_cellchat)

# Calculate the number of links for each object in the list
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})

# Control the dot size in the different datasets
weight.MinMax <- c(min(num.link), max(num.link))

# Generate the signaling role scatter plots for both Het and KO
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

# Combine the plots for easy comparison
combined_plot <- patchwork::wrap_plots(plots = gg)

# Save the combined plot as a PNG file
png(filename = paste0(comparison_path, "signalingRole_scatter_comparison_het_ko.png"), 
    width = 10,
    height = 5, 
    units = "in",
    res = 500)
print(combined_plot)
dev.off()

# ===========================================
# Generate Differential Interaction Plots
# ===========================================
# Merge the two CellChat objects into one for comparison
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))


# Generate the differential interaction plot and show it interactively first
plot_diff <- netVisual_diffInteraction(cellchat_combined, weight.scale = TRUE, title.name = "Difference in Interaction Number (Het vs KO)")

# Save the plot to a PNG file
png(filename = paste0(comparison_path, "diffInteraction_het_ko.png"),
    width = 8,
    height = 10, 
    units = "in",
    res = 300)

# Print the plot to ensure it is saved
print(plot_diff)

# Close the graphics device
dev.off()


# Generate the differential interaction strength plot and show it interactively first
plot_diff_strength <- netVisual_diffInteraction(cellchat_combined, weight.scale = TRUE, measure = "weight", title.name = "Difference in Interaction Strength (Het vs KO)")

# Save the plot to a PNG file
png(filename = paste0(comparison_path, "diffStrength_het_ko.png"),
    width = 8, 
    height = 10,  # Increased height to prevent title cut-off
    units = "in", 
    res = 300)

# Print the plot to ensure it is saved
print(plot_diff_strength)

# Close the graphics device
dev.off()

# ===========================================
# Generate Heatmap Showing Differential Number of Interactions or Interaction Strength
# ===========================================
gg1 <- netVisual_heatmap(cellchat_combined, title.name = "Differential Number of Interactions (Het vs KO)")
gg2 <- netVisual_heatmap(cellchat_combined, measure = "weight", title.name = "Differential Interaction Strength (Het vs KO)")

# Combine the heatmap plots for comparison
combined_heatmap <- gg1 + gg2

# Save the heatmap as a PNG file
png(filename = paste0(comparison_path, "differential_interactions_heatmap_het_ko.png"),
    width = 14, height = 7, units = "in", res = 300)
print(combined_heatmap)
dev.off()



## informationFlow_comparison_het_ko
# Load required packages
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Define paths
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko/")

# Load the saved CellChat objects for Het and KO
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het/het_cellchat.rds"))
ko_cellchat <- readRDS(file = paste0(data_path, "cellchat_ko/ko_cellchat.rds"))

# Merge the two CellChat objects into one for comparison
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# Compare the overall information flow of each signaling pathway or ligand-receptor pair
# Rank the signaling pathways based on the information flow between Het and KO
# This will generate stacked and unstacked bar plots

gg1 <- rankNet(cellchat_combined, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE, font.size = 10) +
  ggtitle("Relative Information Flow: Het vs KO") +
  theme(plot.title = element_text(size = 18, face = "bold"))

gg2 <- rankNet(cellchat_combined, mode = "comparison", measure = "weight", stacked = FALSE, do.stat = TRUE, font.size = 10) +
  ggtitle("Information Flow: Het vs KO") +
  theme(plot.title = element_text(size = 18, face = "bold"))

# Display the plots together for comparison
combined_plot <- gg1 + gg2

# Show the plots interactively before saving
print(combined_plot)

# Save the combined plot as a PNG file
png(filename = paste0(comparison_path, "informationFlow_comparison_het_ko.png"),
    width = 16,
    height = 8, 
    units = "in",
    res = 300)

print(combined_plot)  # Ensure the plot is drawn inside the PNG file
dev.off()



## Identify signaling groups based on their functional similarity between het and ko
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

# Compute network centrality scores for each dataset
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# ===========================================
# Merge the Two CellChat Objects for Comparison
# ===========================================
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# ===========================================
# Identify Signaling Groups Based on Functional Similarity
# ===========================================
ptm <- Sys.time()

# Compute network similarity based on functional similarity
cellchat_combined <- computeNetSimilarityPairwise(cellchat_combined, type = "functional")

# Perform manifold learning of the signaling networks
cellchat_combined <- netEmbedding(cellchat_combined, type = "functional")

# Perform classification learning of the signaling networks
cellchat_combined <- netClustering(cellchat_combined, type = "functional")

# ===========================================
# Visualization in 2D-space (Functional Similarity)
# ===========================================
gg_similarity <- netVisual_embeddingPairwise(cellchat_combined, type = "functional", label.size = 3.5) +
  ggtitle("Functional Similarity between Signaling Networks: Het vs KO") +
  theme(plot.title = element_text(size = 18, face = "bold"))

# Save the visualization as a PNG file
png(filename = paste0(comparison_path, "functional_similarity_het_ko.png"),
    width = 10,
    height = 8, 
    units = "in",
    res = 300)

print(gg_similarity)
dev.off()

# Calculate execution time
execution_time <- Sys.time() - ptm
print(paste("Execution time:", execution_time))

# ===========================================
# Check Data Consistency: Ensure Het and KO Have Sufficient Overlap
# ===========================================
signaling_het <- het_cellchat@netP$pathways
signaling_ko <- ko_cellchat@netP$pathways

intersected_pathways <- intersect(signaling_het, signaling_ko)
print(intersected_pathways)

# ===========================================
# Compute and Visualize Pathway Functional Distance
# ===========================================
distance_similarity <- rankSimilarity(cellchat_combined, type = "functional")

# Add title and increase axis label size
distance_similarity <- distance_similarity +
  ggtitle("Distance Functional Similarity between Signaling Networks: Het vs KO") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 25),  
    axis.title.y = element_text(size = 30)   
  )

print(distance_similarity)

# Save the distance similarity plot
png(filename = paste0(comparison_path, "distance_functional_similarity_het_ko.png"),
    width = 12,
    height = 8, 
    units = "in",
    res = 300)
print(distance_similarity)
dev.off()



## Identify signaling groups based on their structure similarity between het and ko
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

# Compute network centrality scores for each dataset
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# ===========================================
# Merge the Two CellChat Objects for Comparison
# ===========================================
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# ===========================================
# Identify Signaling Groups Based on Structural Similarity
# ===========================================
ptm <- Sys.time()

# Compute network similarity based on structural similarity
cellchat_combined <- computeNetSimilarityPairwise(cellchat_combined, type = "structural")

# Perform manifold learning of the signaling networks based on structural similarity
cellchat_combined <- netEmbedding(cellchat_combined, type = "structural")

# Perform classification learning of the signaling networks
cellchat_combined <- netClustering(cellchat_combined, type = "structural")

# ===========================================
# Visualization in 2D-space (Structural Similarity)
# ===========================================
gg_structural_similarity <- netVisual_embeddingPairwise(cellchat_combined, type = "structural", label.size = 3.5) +
  ggtitle("Structural Similarity between Signaling Networks: Het vs KO") +
  theme(plot.title = element_text(size = 18, face = "bold"))

# Save the visualization as a PNG file
png(filename = paste0(comparison_path, "structural_similarity_het_ko.png"),
    width = 10,
    height = 8, 
    units = "in",
    res = 300)

print(gg_structural_similarity)
dev.off()

# ===========================================
# Visualization with Zoom-in
# ===========================================
gg_structural_similarity_zoom <- netVisual_embeddingPairwiseZoomIn(cellchat_combined, type = "structural", nCol = 2) +
  ggtitle("Zoomed-in Structural Similarity between Signaling Networks: Het vs KO") +
  theme(plot.title = element_text(size = 18, face = "bold"))

# Save the zoomed-in visualization as a PNG file
png(filename = paste0(comparison_path, "structural_similarity_het_ko_zoom.png"),
    width = 12,
    height = 8, 
    units = "in",
    res = 300)

print(gg_structural_similarity_zoom)
dev.off()

# ===========================================
# Compute and Visualize Pathway Structural Distance
# ===========================================
distance_similarity <- rankSimilarity(cellchat_combined, type = "structural")

# Add title and increase axis label size
distance_similarity <- distance_similarity +
  ggtitle("Distance Structural Similarity between Signaling Networks: Het vs KO") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 25),  
    axis.title.y = element_text(size = 30)   
  )

print(distance_similarity)

# Save the distance similarity plot
png(filename = paste0(comparison_path, "distance_structural_similarity_het_ko.png"),
    width = 12,
    height = 8, 
    units = "in",
    res = 300)
print(distance_similarity)
dev.off()

# Calculate execution time
execution_time <- Sys.time() - ptm
print(paste("Execution time:", execution_time))




## heatmap_outgoing_comparison_het_ko
# ===========================================
# Load Required Libraries
# ===========================================
library(CellChat)
library(ComplexHeatmap)

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

# Compute centrality scores for each dataset
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# ===========================================
# Merge the Two CellChat Objects for Comparison
# ===========================================
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)

# ===========================================
# Generate Heatmaps for Comparing Incoming Signaling Patterns
# ===========================================
ht1_incoming <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[1], width = 8, height = 15, font.size.title = 14)
ht2_incoming <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[2], width = 8, height = 15, font.size.title = 14)

# Save the combined incoming signaling heatmap
png(filename = paste0(comparison_path, "heatmap_incoming_comparison_het_ko.png"),
    width = 20,
    height = 20,
    units = "in",
    res = 500)

draw(ht1_incoming + ht2_incoming, ht_gap = unit(0.5, "cm"))
dev.off()

# ===========================================
# Generate Heatmaps for Comparing Outgoing Signaling Patterns
# ===========================================
ht1_outgoing <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[1], width = 8, height = 15, font.size.title = 14)
ht2_outgoing <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[2], width = 8, height = 15, font.size.title = 14)

# Save the combined outgoing signaling heatmap
png(filename = paste0(comparison_path, "heatmap_outgoing_comparison_het_ko.png"),
    width = 20,
    height = 20,
    units = "in",
    res = 500)

draw(ht1_outgoing + ht2_outgoing, ht_gap = unit(0.5, "cm"))
dev.off()




# ============================================================================
# Compare the Total Number of Interactions and Interaction Strength (Het vs KO)
# ============================================================================
# This script compares the total number of cell-cell interactions and 
# the interaction strength between Het and KO samples using CellChat. 
# It visualizes whether cell-cell communication is enhanced or reduced 
# in different conditions.
# ============================================================================

# Load required packages
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggrepel)

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

# Compute centrality scores for each dataset
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# ===========================================
# Merge the Two CellChat Objects for Comparison
# ===========================================
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# ===========================================
# Compare the Total Number of Interactions
# ===========================================
gg1 <- compareInteractions(cellchat_combined, show.legend = FALSE, group = c(1,2), size.text = 15, color.use = c("#66B2FF", "#E63E62"))

# Compare the Interaction Strength
gg2 <- compareInteractions(cellchat_combined, show.legend = FALSE, group = c(1,2), size.text = 15, measure = "weight", color.use = c("#66B2FF", "#E63E62"))

# Combine the plots and move the title higher using `margin`
combined_plot <- gg1 + gg2 + 
  plot_annotation(title = "Compare the Total Number of Interactions and Interaction Strength (Het vs KO)") & 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", margin = margin(b = 40)))  # Increase `b` value for more space

# Print the plot
print(combined_plot)

# Save the plot as a PNG file
png(filename = paste0(comparison_path, "compare_interactions_het_ko.png"), 
    width = 5, height = 5, units = "in", res = 500)

print(combined_plot)
dev.off()





# ==============================================================================
# Generate Violin Plots for Ligand-Receptor Pairs Comparing Het vs KO
# ==============================================================================

# Load required packages
library(CellChat)
library(ggplot2)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko/")
violin_path <- paste0(comparison_path, "violin_LR_comparison/")

# ===========================================
# Load the Saved CellChat Objects
# ===========================================
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het/het_cellchat.rds"))
ko_cellchat <- readRDS(file = paste0(data_path, "cellchat_ko/ko_cellchat.rds"))

# ===========================================
# Merge the Two CellChat Objects for Comparison
# ===========================================
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# Add metadata column for dataset labels (Het vs KO)
cellchat_combined@meta$datasets <- factor(cellchat_combined@meta$datasets, levels = c("Het", "KO"))

# ===========================================
# Identify Common Pathways in Both Conditions
# ===========================================
pathways.show.all <- intersect(het_cellchat@netP$pathways, ko_cellchat@netP$pathways)

# ===========================================
# Generate Violin Plots for Ligand-Receptor Pairs
# ===========================================
for (i in 1:length(pathways.show.all)) {
  pathway <- pathways.show.all[i]
  
  # Generate violin plot
  gg <- plotGeneExpression(cellchat_combined, 
                           signaling = pathway, 
                           enriched.only = TRUE, 
                           type = "violin", 
                           split.by = "datasets",  # Differentiating Het and KO
                           colors.ggplot = TRUE) +
        ggtitle(paste("Expression Levels for Signaling Pathway Components:", pathway))

  # Save the violin plot
  ggsave(filename = paste0(violin_path, "violin_LR_comparison_", pathway, "_het_ko.png"),
         plot = gg, width = 12, height = 7, dpi = 300)
}




#### remove Red_Blood_Cells
# ===========================================
# Load Required Libraries
# ===========================================
rm(list = ls())

library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(dplyr)

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

# Check Data Consistency: Ensure both have sufficient overlap in signaling pathways and cell types
signaling_het <- het_cellchat@netP$pathways
signaling_ko <- ko_cellchat@netP$pathways

# Prepare the CellChat objects by aggregating the communication networks
het_cellchat <- aggregateNet(het_cellchat)
ko_cellchat <- aggregateNet(ko_cellchat)

# Calculate the centrality scores for the signaling roles
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# Create a list of CellChat objects for comparison
object.list <- list(Het = het_cellchat, KO = ko_cellchat)

# Calculate the number of links for each object in the list
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})

# Control the dot size in the different datasets
weight.MinMax <- c(min(num.link), max(num.link))

# Generate the signaling role scatter plots for both Het and KO
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, font.size = 14, font.size.title = 14)
}

# Combine the plots for easy comparison
combined_plot <- patchwork::wrap_plots(plots = gg)

# Save the combined plot as a PNG file
png(filename = paste0(comparison_path, "signalingRole_scatter_comparison_het_ko.png"), 
    width = 10,
    height = 5, 
    units = "in",
    res = 500)
print(combined_plot)
dev.off()



# ===========================================
# Scatter Plot of Incoming Interaction Strengths:
# Comparing Het vs KO CellChat Results
## absoult_value_scatter_incoming_het_vs_ko
# ===========================================
# ===========================================
# Load Required Libraries
# ===========================================
# rm(list = ls())  # Clear environment
library(CellChat)       # Cell-cell communication analysis
library(ggplot2)        # Plotting
library(ggrepel)        # Better text labels in ggplot
library(dplyr)          # Data manipulation

# ===========================================
# Define File Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko_02/")

# Create output directory if it doesn't exist
if (!dir.exists(comparison_path)) dir.create(comparison_path, recursive = TRUE)

# ===========================================
# Load CellChat Objects for Het and KO Conditions
# ===========================================
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het_02/het_cellchat_02.rds"))
ko_cellchat  <- readRDS(file = paste0(data_path, "cellchat_ko_02/ko_cellchat_02.rds"))

# ===========================================
# Preprocess: Aggregate Network & Compute Centrality Measures
# ===========================================
het_cellchat <- aggregateNet(het_cellchat)  # Summarize communication networks
ko_cellchat  <- aggregateNet(ko_cellchat)

## scatter_incoming_het_vs_ko.png
# ===========================================
# Compute Incoming Interaction Strengths
# (Aggregate over all pathways)
# ===========================================
incoming_het <- rowSums(het_cellchat@netP$prob, dims = 1)
incoming_ko  <- rowSums(ko_cellchat@netP$prob, dims = 1)


# Make sure the dimensions match
common_celltypes <- intersect(names(incoming_het), names(incoming_ko))

df_compare <- data.frame(
  CellType = common_celltypes,
  incoming_Het = incoming_het[common_celltypes],
  incoming_KO  = incoming_ko[common_celltypes]
)


p <- ggplot(df_compare, aes(x = incoming_Het, y = incoming_KO)) +
  geom_point(size = 4, color = "steelblue") +
  geom_text_repel(aes(label = CellType), size = 4) +
  labs(
    x = "Het: Incoming interaction strength",
    y = "KO: Incoming interaction strength",
    title = "Comparison of Incoming Interaction Strengths: Het vs KO"
  ) +
  theme_minimal(base_size = 14)

# Save the plot
ggsave(
  filename = paste0(comparison_path, "scatter_incoming_het_vs_ko.png"),
  plot = p, width = 6, height = 6, dpi = 500
)

ggsave(paste0(comparison_path, "scatter_incoming_het_vs_ko.png"), plot = p, width = 6, height = 6, dpi = 500)




## bar_incoming_centrality_het_vs_ko_normalized
library(reshape2)

# Compute centrality for both
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat  <- netAnalysis_computeCentrality(ko_cellchat,  slot.name = "netP")

# Extract and average indegree centrality
get_avg_indegree <- function(cellchat_obj) {
  centr_list <- lapply(cellchat_obj@netP$centr, function(x) x$indeg)
  Reduce("+", centr_list) / length(centr_list)
}

incoming_het <- get_avg_indegree(het_cellchat)
incoming_ko  <- get_avg_indegree(ko_cellchat)

# Combine into one data frame
df_compare <- data.frame(
  CellType = names(incoming_het),
  Het = incoming_het,
  KO  = incoming_ko
)

# Convert to long format for ggplot2
df_long <- melt(df_compare, id.vars = "CellType", variable.name = "Condition", value.name = "Centrality")

# Plot
p <- ggplot(df_long, aes(x = reorder(CellType, Centrality), y = Centrality, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(
    title = "Normalized Incoming Centrality: Het vs KO",
    x = "Cell Type",
    y = "Average Indegree Centrality"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("Het" = "steelblue", "KO" = "darkorange"))

# Show plot
print(p)

# Save
ggsave(
  filename = paste0(comparison_path, "bar_incoming_centrality_het_vs_ko_normalized.png"),
  plot = p, width = 8, height = 6, dpi = 500
)


## scatter_incoming_centrality_het_vs_ko
# Compute centrality
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat  <- netAnalysis_computeCentrality(ko_cellchat,  slot.name = "netP")

# Get average incoming (indegree) centrality
get_avg_indegree <- function(cellchat_obj) {
  centr_list <- lapply(cellchat_obj@netP$centr, function(x) x$indeg)
  Reduce("+", centr_list) / length(centr_list)
}

incoming_het <- get_avg_indegree(het_cellchat)
incoming_ko  <- get_avg_indegree(ko_cellchat)

# Ensure matching cell types
common_cells <- intersect(names(incoming_het), names(incoming_ko))
df_scatter <- data.frame(
  CellType = common_cells,
  Het = incoming_het[common_cells],
  KO  = incoming_ko[common_cells]
)

# Plot scatter
p <- ggplot(df_scatter, aes(x = Het, y = KO)) +
  geom_point(color = "steelblue", size = 4) +
  geom_text_repel(aes(label = CellType), size = 4, max.overlaps = 100) +
  labs(
    title = "Scatter Plot of Normalized Incoming Centrality",
    x = "Het: Average Indegree Centrality",
    y = "KO: Average Indegree Centrality"
  ) +
  theme_minimal(base_size = 14)

# Save plot
ggsave(
  filename = paste0(comparison_path, "scatter_incoming_centrality_het_vs_ko.png"),
  plot = p, width = 6, height = 6, dpi = 500
)



## scatter_incoming_centrality_het_vs_ko_colored
p <- ggplot(df_scatter, aes(x = Het, y = KO, color = CellType)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = CellType), size = 4, max.overlaps = 100, show.legend = FALSE) +
  labs(
    title = "Scatter Plot of Normalized Incoming Centrality",
    x = "Het: Average Indegree Centrality",
    y = "KO: Average Indegree Centrality",
    color = "Cell Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",           # Hide legend
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

# Save the plot
ggsave(
  filename = paste0(comparison_path, "scatter_incoming_centrality_het_vs_ko_colored.png"),
  plot = p, width = 7, height = 6, dpi = 500
)


## scatter_incoming_centrality_with_count.png
# Compute number of incoming links per cell type (sum over senders and pathways)
count_het <- apply(het_cellchat@net$count, 2, sum)
count_ko  <- apply(ko_cellchat@net$count, 2, sum)

# Match cell types
common_cells <- intersect(names(count_het), names(count_ko))

# Merge with scatter dataframe
df_scatter$Count <- (count_het[common_cells] + count_ko[common_cells]) / 2  # average count

# Scale x and y to percentage
df_scatter$Het <- df_scatter$Het * 100
df_scatter$KO  <- df_scatter$KO * 100

# Load RColorBrewer for color palettes
library(colorspace)

# Define distinct colors for each unique cell type
celltypes <- unique(df_scatter$CellType)
custom_colors <- setNames(qualitative_hcl(length(celltypes), palette = "Dark 3"), celltypes)

p <- ggplot(df_scatter, aes(x = Het, y = KO, color = CellType, size = Count)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = CellType), size = 4, max.overlaps = 100, box.padding = 0.6, point.padding = 0.6, show.legend = FALSE) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Incoming Interaction Strength",
    x = "Het: Incoming Interaction Strength (%)",
    y = "KO: Incoming Interaction Strength (%)"
  ) +
  scale_size_continuous(range = c(3, 10)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  ) +
  guides(color = "none")  # hide color legend, keep size legend


# Save
ggsave(
  filename = paste0(comparison_path, "scatter_incoming_centrality_with_count.png"),
  plot = p, width = 7, height = 6, dpi = 500
)






# ===========================================
# Generate Differential Interaction Plots
# ===========================================
# Merge the two CellChat objects into one for comparison
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))


# Generate the differential interaction plot and show it interactively first
plot_diff <- netVisual_diffInteraction(cellchat_combined, weight.scale = TRUE, title.name = "Difference in Interaction Number (Het vs KO)", vertex.label.cex = 1.3)

# Save the plot to a PNG file
png(filename = paste0(comparison_path, "diffInteraction_het_ko.png"),
    width = 8,
    height = 10, 
    units = "in",
    res = 300)

# Print the plot to ensure it is saved
print(plot_diff)

# Close the graphics device
dev.off()


# Generate the differential interaction strength plot and show it interactively first
plot_diff_strength <- netVisual_diffInteraction(cellchat_combined, weight.scale = TRUE, measure = "weight", title.name = "Difference in Interaction Strength (Het vs KO)", vertex.label.cex = 1.3)

# Save the plot to a PNG file
png(filename = paste0(comparison_path, "diffStrength_het_ko.png"),
    width = 8, 
    height = 10,  # Increased height to prevent title cut-off
    units = "in", 
    res = 300)

# Print the plot to ensure it is saved
print(plot_diff_strength)

# Close the graphics device
dev.off()

# ===========================================
# Generate Heatmap Showing Differential Number of Interactions or Interaction Strength
# ===========================================
gg1 <- netVisual_heatmap(cellchat_combined, title.name = "Differential Number of Interactions (Het vs KO)", font.size = 15, font.size.title = 15)
gg2 <- netVisual_heatmap(cellchat_combined, measure = "weight", title.name = "Differential Interaction Strength (Het vs KO)", font.size = 15, font.size.title = 15)

# Combine the heatmap plots for comparison
combined_heatmap <- gg1 + gg2

# Save the heatmap as a PNG file
png(filename = paste0(comparison_path, "differential_interactions_heatmap_het_ko.png"),
    width = 14, height = 7, units = "in", res = 300)
print(combined_heatmap)
dev.off()



## informationFlow_comparison_het_ko
# Load required packages
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Define paths
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko_02/")

# Load the saved CellChat objects for Het and KO
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het_02/het_cellchat_02.rds"))
ko_cellchat <- readRDS(file = paste0(data_path, "cellchat_ko_02/ko_cellchat_02.rds"))

# Merge the two CellChat objects into one for comparison
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# Compare the overall information flow of each signaling pathway or ligand-receptor pair
# Rank the signaling pathways based on the information flow between Het and KO
# This will generate stacked and unstacked bar plots

gg1 <- rankNet(cellchat_combined, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE, font.size = 10) +
  ggtitle("Relative Information Flow: Het vs KO") +
  theme(plot.title = element_text(size = 18, face = "bold"))

gg2 <- rankNet(cellchat_combined, mode = "comparison", measure = "weight", stacked = FALSE, do.stat = TRUE, font.size = 10) +
  ggtitle("Information Flow: Het vs KO") +
  theme(plot.title = element_text(size = 18, face = "bold"))

# Display the plots together for comparison
combined_plot <- gg1 + gg2

# Show the plots interactively before saving
print(combined_plot)

# Save the combined plot as a PNG file
png(filename = paste0(comparison_path, "informationFlow_comparison_het_ko.png"),
    width = 16,
    height = 8, 
    units = "in",
    res = 300)

print(combined_plot)  # Ensure the plot is drawn inside the PNG file
dev.off()





## Compare the overall information flow of top structurally different signaling pathways (Het vs KO)
# Load required packages
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Define paths
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko_02/")

# Load the saved CellChat objects for Het and KO
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het_02/het_cellchat_02.rds"))
ko_cellchat <- readRDS(file = paste0(data_path, "cellchat_ko_02/ko_cellchat_02.rds"))

# Merge the two CellChat objects into one for comparison
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# Define pathways based on structural distance difference from the image
selected_signaling <- c("GAP", "CNTN", "PSAP", "PCDH", "NCAM", "CADM", 
                        "CypA", "PECAM1", "FLRT", "ApoE", "PTN", "APP", 
                        "ESAM", "CDH5", "EDN", "MK")

# Generate stacked and unstacked bar plots for these selected pathways
gg1 <- rankNet(cellchat_combined, mode = "comparison", measure = "weight", 
               stacked = TRUE, do.stat = TRUE, font.size = 10, 
               signaling = selected_signaling, color.use = c("#66B2FF", "#E63E62")) +
  ggtitle("Relative Information Flow: Top Structurally Different Pathways (Het vs KO)") +
  theme(plot.title = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(color = "black"))

gg2 <- rankNet(cellchat_combined, mode = "comparison", measure = "weight", 
               stacked = FALSE, do.stat = TRUE, font.size = 10, 
               signaling = selected_signaling, color.use = c("#66B2FF", "#E63E62")) +
  ggtitle("Information Flow: Top Structurally Different Pathways (Het vs KO)") +
  theme(plot.title = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(color = "black"))

# Combine and display the plots
combined_plot <- gg1 + gg2
print(combined_plot)

# Save the plot as a PNG
png(filename = paste0(comparison_path, "informationFlow_structuralPathways_het_ko.png"),
    width = 16,
    height = 8,
    units = "in",
    res = 300)
print(combined_plot)
dev.off()


gg2 <- rankNet(cellchat_combined, 
               mode = "comparison", 
               measure = "weight", 
               stacked = FALSE, 
               do.stat = TRUE, 
               signaling = selected_signaling, 
               bar.w = 0.98, 
               do.flip = FALSE, 
               color.use = c("#66B2FF", "#E63E62")) +
  ggtitle("Information Flow: Structurally Different Pathways (Het vs KO)") +
  theme(
    plot.title = element_text(size = 8, face = "bold"),
    axis.title.y = element_text(size = 12),     # Y-axis label
    axis.text.y = element_text(size = 12, color = "black"),    # Y-axis tick labels
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, color = "black"), # X-axis tick labels
    legend.title = element_text(size = 10),     # Legend title
    legend.text = element_text(size = 8),                     # Legend entries
    legend.position = c(0.72, 0.85),  # Adjust horizontal to be over GAP
    legend.justification = c("left", "top")
  )

# Save only gg2 (unstacked bar plot)
png(filename = paste0(comparison_path, "informationFlow_structuralPathways_gg2_het_ko.png"),
    width = 4,
    height = 3,
    units = "in",
    res = 300)

print(gg2)  # Only save the gg2 plot
dev.off()



## Identify signaling groups based on their functional similarity between het and ko
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

# Compute network centrality scores for each dataset
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# ===========================================
# Merge the Two CellChat Objects for Comparison
# ===========================================
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# ===========================================
# Identify Signaling Groups Based on Functional Similarity
# ===========================================
ptm <- Sys.time()

# Set random seed for reproducibility
set.seed(2025)


# Compute network similarity based on functional similarity
cellchat_combined <- computeNetSimilarityPairwise(cellchat_combined, type = "functional")

# Perform manifold learning of the signaling networks
cellchat_combined <- netEmbedding(cellchat_combined, type = "functional")

# Perform classification learning of the signaling networks
cellchat_combined <- netClustering(cellchat_combined, type = "functional")

# ===========================================
# Visualization in 2D-space (Functional Similarity)
# ===========================================
gg_similarity <- netVisual_embeddingPairwise(cellchat_combined, type = "functional", label.size = 3.5) +
  ggtitle("Functional Similarity between Signaling Networks: Het vs KO") +
  theme(plot.title = element_text(size = 18, face = "bold"))

# Save the visualization as a PNG file
png(filename = paste0(comparison_path, "functional_similarity_het_ko.png"),
    width = 10,
    height = 8, 
    units = "in",
    res = 300)

print(gg_similarity)
dev.off()

# Calculate execution time
execution_time <- Sys.time() - ptm
print(paste("Execution time:", execution_time))

# ===========================================
# Check Data Consistency: Ensure Het and KO Have Sufficient Overlap
# ===========================================
signaling_het <- het_cellchat@netP$pathways
signaling_ko <- ko_cellchat@netP$pathways

intersected_pathways <- intersect(signaling_het, signaling_ko)
print(intersected_pathways)

# ===========================================
# Compute and Visualize Pathway Functional Distance
# ===========================================
distance_similarity <- rankSimilarity(cellchat_combined, type = "functional")

# Add title and increase axis label size
distance_similarity <- distance_similarity +
  ggtitle("Distance Functional Similarity between Signaling Networks: Het vs KO") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 25),  
    axis.title.y = element_text(size = 30)   
  )

print(distance_similarity)

# Save the distance similarity plot
png(filename = paste0(comparison_path, "distance_functional_similarity_het_ko_04.png"),
    width = 12,
    height = 8, 
    units = "in",
    res = 300)
print(distance_similarity)
dev.off()



## Identify signaling groups based on their structure similarity between het and ko
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

# Compute network centrality scores for each dataset
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# ===========================================
# Merge the Two CellChat Objects for Comparison
# ===========================================
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# ===========================================
# Identify Signaling Groups Based on Structural Similarity
# ===========================================
ptm <- Sys.time()

# Set random seed for reproducibility
set.seed(2027)

# Compute network similarity based on structural similarity
cellchat_combined <- computeNetSimilarityPairwise(cellchat_combined, type = "structural")

# Perform manifold learning of the signaling networks based on structural similarity
cellchat_combined <- netEmbedding(cellchat_combined, type = "structural")

# Perform classification learning of the signaling networks
cellchat_combined <- netClustering(cellchat_combined, type = "structural")

# ===========================================
# Visualization in 2D-space (Structural Similarity)
# ===========================================
gg_structural_similarity <- netVisual_embeddingPairwise(cellchat_combined, type = "structural", label.size = 3.5) +
  ggtitle("Structural Similarity between Signaling Networks: Het vs KO") +
  theme(plot.title = element_text(size = 18, face = "bold"))

# Save the visualization as a PNG file
png(filename = paste0(comparison_path, "structural_similarity_het_ko.png"),
    width = 10,
    height = 8, 
    units = "in",
    res = 300)

print(gg_structural_similarity)
dev.off()

# ===========================================
# Visualization with Zoom-in
# ===========================================
gg_structural_similarity_zoom <- netVisual_embeddingPairwiseZoomIn(cellchat_combined, type = "structural", nCol = 2) +
  ggtitle("Zoomed-in Structural Similarity between Signaling Networks: Het vs KO") +
  theme(plot.title = element_text(size = 18, face = "bold"))

# Save the zoomed-in visualization as a PNG file
png(filename = paste0(comparison_path, "structural_similarity_het_ko_zoom.png"),
    width = 12,
    height = 8, 
    units = "in",
    res = 300)

print(gg_structural_similarity_zoom)
dev.off()

# ===========================================
# Compute and Visualize Pathway Structural Distance
# ===========================================
distance_similarity <- rankSimilarity(cellchat_combined, type = "structural", font.size = 15)

# Add title and increase axis label size
distance_similarity <- distance_similarity +
  ggtitle("Distance Structural Similarity between Signaling Networks: Het vs KO") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 25),  
    axis.title.y = element_text(size = 30)   
  )

print(distance_similarity)

# Save the distance similarity plot
png(filename = paste0(comparison_path, "distance_structural_similarity_het_ko_02.png"),
    width = 12,
    height = 8, 
    units = "in",
    res = 300)
print(distance_similarity)
dev.off()

# Calculate execution time
execution_time <- Sys.time() - ptm
print(paste("Execution time:", execution_time))





# ===========================================
# Identify signaling groups based on their structure similarity between Het and KO
# Extract data from ggplot object
# ===========================================

# Load libraries
library(CellChat)
library(ggplot2)
library(dplyr)

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

# Compute network centrality scores for each dataset
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# ===========================================
# Merge the Two CellChat Objects for Comparison
# ===========================================
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# ===========================================
# Identify Signaling Groups Based on Structural Similarity
# ===========================================

# Set random seed for reproducibility
#set.seed(12345)
set.seed(12346)



# Compute network similarity based on structural similarity
cellchat_combined <- computeNetSimilarityPairwise(cellchat_combined, type = "structural")

# Manifold learning and clustering of signaling networks
cellchat_combined <- netEmbedding(cellchat_combined, type = "structural")
cellchat_combined <- netClustering(cellchat_combined, type = "structural")

# ===========================================
# Compute and Visualize Pathway Structural Distance
# ===========================================
distance_similarity <- rankSimilarity(cellchat_combined, type = "structural")

# Extract data from ggplot object
distance_similarity_data <- distance_similarity$data

# ===========================================
# Create Bar Plot for All Pathways
# ===========================================
distance_similarity_all <- ggplot(distance_similarity_data, aes(y = dist, x = reorder(name, -dist))) +
  geom_bar(stat = "identity", fill = "#4B9CD3") +  # Custom blue
  labs(y = "Pathway Distance", x = "Signaling Pathway") + 
  ggtitle("Distance Structural Similarity :\nHet vs KO") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    aspect.ratio = 1.5
  )

# Show the plot
print(distance_similarity_all)

# ===========================================
# Save the Plot to File
# ===========================================
png(filename = paste0(comparison_path, "distance_structural_similarity_het_ko_06.png"),
    width = 10,
    height = 6,
    units = "in",
    res = 300)

print(distance_similarity_all)
dev.off()






## heatmap_outgoing_comparison_het_ko
# ===========================================
# Load Required Libraries
# ===========================================
library(CellChat)
library(ComplexHeatmap)

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

# Compute centrality scores for each dataset
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# ===========================================
# Merge the Two CellChat Objects for Comparison
# ===========================================
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)

# ===========================================
# Generate Heatmaps for Comparing Incoming Signaling Patterns
# ===========================================
ht1_incoming <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[1], width = 8, height = 15, font.size.title = 14, color.heatmap = "Reds")
ht2_incoming <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[2], width = 8, height = 15, font.size.title = 14, color.heatmap = "Reds")

# Save the combined incoming signaling heatmap
png(filename = paste0(comparison_path, "heatmap_incoming_comparison_het_ko_02.png"),
    width = 20,
    height = 20,
    units = "in",
    res = 500)

draw(ht1_incoming + ht2_incoming, ht_gap = unit(0.5, "cm"))
dev.off()

# ===========================================
# Generate Heatmaps for Comparing Outgoing Signaling Patterns
# ===========================================
ht1_outgoing <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[1], width = 8, height = 15, font.size.title = 14, color.heatmap = "Reds")
ht2_outgoing <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[2], width = 8, height = 15, font.size.title = 14, color.heatmap = "Reds")

# Save the combined outgoing signaling heatmap
png(filename = paste0(comparison_path, "heatmap_outgoing_comparison_het_ko_02.png"),
    width = 20,
    height = 20,
    units = "in",
    res = 500)

draw(ht1_outgoing + ht2_outgoing, ht_gap = unit(0.5, "cm"))
dev.off()



## heatmap_incoming_outgoing_comparison_het_ko selected pathways
# ===========================================
# Load Required Libraries
# ===========================================
library(CellChat)
library(ComplexHeatmap)
library(grid)  # For 'unit' used in draw()

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko_02/")

# ===========================================
# Load the Saved CellChat Objects
# ===========================================
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het_02/het_cellchat_02.rds"))
ko_cellchat  <- readRDS(file = paste0(data_path, "cellchat_ko_02/ko_cellchat_02.rds"))

# ===========================================
# Aggregate Communication Networks & Compute Centrality
# ===========================================
het_cellchat <- aggregateNet(het_cellchat)
ko_cellchat  <- aggregateNet(ko_cellchat)

het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat  <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# ===========================================
# Merge CellChat Objects for Comparison
# ===========================================
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)

# ===========================================
# Define Selected Pathways (based on figure)
# ===========================================
selected.pathways <- c(
  "GAP", "CNTN", "PSAP", "PCDH", "NCAM", "CADM", 
  "CypA", "PECAM1", "FLRT", "ApoE", "PTN", 
  "APP", "ESAM", "CDH5", "EDN", "MK"
)

# Keep only pathways available in the data
selected.pathways <- intersect(selected.pathways, pathway.union)

# ===========================================
# Generate Incoming Heatmaps
# ===========================================
ht1_incoming <- netAnalysis_signalingRole_heatmap(
  object.list[[1]], pattern = "incoming", signaling = selected.pathways,
  title = names(object.list)[1], width = 8, height = 10, 
  font.size = 14, font.size.title = 14, color.heatmap = "Reds"
)

ht2_incoming <- netAnalysis_signalingRole_heatmap(
  object.list[[2]], pattern = "incoming", signaling = selected.pathways,
  title = names(object.list)[2], width = 8, height = 10, 
  font.size = 14, font.size.title = 14, color.heatmap = "Reds"
)

# Save Incoming Heatmap
png(filename = paste0(comparison_path, "heatmap_incoming_comparison_selected_pathways_het_ko.png"),
    width = 24, height = 20, units = "in", res = 500)
draw(ht1_incoming + ht2_incoming, ht_gap = unit(0.5, "cm"))
dev.off()



## heatmap_incoming_outgoing_comparison_het_ko selected pathways_02
# ===========================================
# Load Required Libraries
# ===========================================
library(CellChat)
library(ComplexHeatmap)
library(grid)  # For 'unit' used in draw()

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko_02/")

# ===========================================
# Load the Saved CellChat Objects
# ===========================================
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het_02/het_cellchat_02.rds"))
ko_cellchat  <- readRDS(file = paste0(data_path, "cellchat_ko_02/ko_cellchat_02.rds"))

# ===========================================
# Aggregate Communication Networks & Compute Centrality
# ===========================================
het_cellchat <- aggregateNet(het_cellchat)
ko_cellchat  <- aggregateNet(ko_cellchat)

het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat  <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# ===========================================
# Merge CellChat Objects for Comparison
# ===========================================
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)

# ===========================================
# Define Selected Pathways (based on figure)
# ===========================================
selected.pathways <- c(
  "APP", "PCDH", "CADM" 
  )

# Keep only pathways available in the data
selected.pathways <- intersect(selected.pathways, pathway.union)

# ===========================================
# Generate Incoming Heatmaps
# ===========================================
ht1_incoming <- netAnalysis_signalingRole_heatmap(
  object.list[[1]], pattern = "incoming", signaling = selected.pathways,
  title = names(object.list)[1], width = 8, height = 2, 
  font.size = 14, font.size.title = 14, color.heatmap = "Reds"
)

ht2_incoming <- netAnalysis_signalingRole_heatmap(
  object.list[[2]], pattern = "incoming", signaling = selected.pathways,
  title = names(object.list)[2], width = 8, height = 2, 
  font.size = 14, font.size.title = 14, color.heatmap = "Reds"
)

# Save Incoming Heatmap
png(filename = paste0(comparison_path, "heatmap_incoming_comparison_selected_pathways_het_ko_02.png"),
    width = 24, height = 20, units = "in", res = 500)
draw(ht1_incoming + ht2_incoming, ht_gap = unit(0.5, "cm"))
dev.off()


# ===========================================
# Generate Outgoing Heatmaps
# ===========================================
ht1_outgoing <- netAnalysis_signalingRole_heatmap(
  object.list[[1]], pattern = "outgoing", signaling = selected.pathways,
  title = names(object.list)[1], width = 8, height = 2, 
  font.size = 14, font.size.title = 14, color.heatmap = "Reds"
)

ht2_outgoing <- netAnalysis_signalingRole_heatmap(
  object.list[[2]], pattern = "outgoing", signaling = selected.pathways,
  title = names(object.list)[2], width = 8, height = 2, 
  font.size = 14, font.size.title = 14, color.heatmap = "Reds"
)

# Save Outgoing Heatmap
png(filename = paste0(comparison_path, "heatmap_outgoing_comparison_selected_pathways_het_ko_02.png"),
    width = 24, height = 20, units = "in", res = 500)
draw(ht1_outgoing + ht2_outgoing, ht_gap = unit(0.5, "cm"))
dev.off()











# ============================================================================
# Compare the Total Number of Interactions and Interaction Strength (Het vs KO)
# ============================================================================
# This script compares the total number of cell-cell interactions and 
# the interaction strength between Het and KO samples using CellChat. 
# It visualizes whether cell-cell communication is enhanced or reduced 
# in different conditions.
# ============================================================================

# Load required packages
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggrepel)

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

# Compute centrality scores for each dataset
het_cellchat <- netAnalysis_computeCentrality(het_cellchat, slot.name = "netP")
ko_cellchat <- netAnalysis_computeCentrality(ko_cellchat, slot.name = "netP")

# ===========================================
# Merge the Two CellChat Objects for Comparison
# ===========================================
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# ===========================================
# Compare the Total Number of Interactions
# ===========================================
gg1 <- compareInteractions(cellchat_combined, show.legend = FALSE, group = c(1,2), size.text = 15, color.use = c("#66B2FF", "#E63E62"))

# Compare the Interaction Strength
gg2 <- compareInteractions(cellchat_combined, show.legend = FALSE, group = c(1,2), size.text = 15, measure = "weight", color.use = c("#66B2FF", "#E63E62"))

# Combine the plots and move the title higher using `margin`
combined_plot <- gg1 + gg2 + 
  plot_annotation(title = "Compare the Total Number of Interactions and Interaction Strength (Het vs KO)") & 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", margin = margin(b = 40)))  # Increase `b` value for more space

# Print the plot
print(combined_plot)

# Save the plot as a PNG file
png(filename = paste0(comparison_path, "compare_interactions_het_ko.png"), 
    width = 5, height = 5, units = "in", res = 500)

print(combined_plot)
dev.off()





# ==============================================================================
# Generate Violin Plots for Ligand-Receptor Pairs Comparing Het vs KO
# ==============================================================================

# Load required packages
library(CellChat)
library(ggplot2)

# ===========================================
# Define Paths
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
comparison_path <- paste0(data_path, "comparasion_het_ko_02/")
violin_path <- paste0(comparison_path, "violin_LR_comparison/")

# ===========================================
# Load the Saved CellChat Objects
# ===========================================
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het_02/het_cellchat_02.rds"))
ko_cellchat <- readRDS(file = paste0(data_path, "cellchat_ko_02/ko_cellchat_02.rds"))

# ===========================================
# Merge the Two CellChat Objects for Comparison
# ===========================================
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# Add metadata column for dataset labels (Het vs KO)
cellchat_combined@meta$datasets <- factor(cellchat_combined@meta$datasets, levels = c("Het", "KO"))

# ===========================================
# Identify Common Pathways in Both Conditions
# ===========================================
pathways.show.all <- intersect(het_cellchat@netP$pathways, ko_cellchat@netP$pathways)

# ===========================================
# Generate Violin Plots for Ligand-Receptor Pairs
# ===========================================
for (i in 1:length(pathways.show.all)) {
  pathway <- pathways.show.all[i]
  
  # Generate violin plot
  gg <- plotGeneExpression(cellchat_combined, 
                           signaling = pathway, 
                           enriched.only = TRUE, 
                           type = "violin", 
                           split.by = "datasets",  # Differentiating Het and KO
                           colors.ggplot = TRUE) +
        ggtitle(paste("Expression Levels for Signaling Pathway Components:", pathway))

  # Save the violin plot
  ggsave(filename = paste0(violin_path, "violin_LR_comparison_", pathway, "_het_ko.png"),
         plot = gg, width = 12, height = 7, dpi = 300)
}





# ====================================================
# Wilcoxon p-value testing for selected gene-cell types (no pathways, no plots)
# ====================================================

library(CellChat)
library(dplyr)

# Define paths
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
output_path <- paste0(data_path, "comparasion_het_ko_02/")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# Load and merge CellChat objects
het_cellchat <- readRDS(file = paste0(data_path, "cellchat_het_02/het_cellchat_02.rds"))
ko_cellchat <- readRDS(file = paste0(data_path, "cellchat_ko_02/ko_cellchat_02.rds"))
object.list <- list(Het = het_cellchat, KO = ko_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))
cellchat_combined@meta$datasets <- factor(cellchat_combined@meta$datasets, levels = c("Het", "KO"))

# Define gene–cell type combinations
tests_to_run <- list(
  list(gene = "App",     cell_type = "Stem_A"),
  list(gene = "App",     cell_type = "Stem_B"),
  list(gene = "Pcdhga9", cell_type = "Stem_A"),
  list(gene = "Pcdhga9", cell_type = "Stem_B"),
  list(gene = "Cadm1",   cell_type = "Stem_A"),
  list(gene = "Cadm1",   cell_type = "Stem_B"),
  list(gene = "Cadm1",   cell_type = "Mitotic_B"),
  list(gene = "Cadm1",   cell_type = "Oligodendrocytes"),
  list(gene = "Cadm1",   cell_type = "Purkinje_Cells")
)

# Prepare results table
pval_results <- data.frame(
  Gene = character(),
  CellType = character(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)

# Run tests
expr_matrix <- cellchat_combined@data.signaling
meta <- cellchat_combined@meta

for (test_case in tests_to_run) {
  gene <- test_case$gene
  cell_type <- test_case$cell_type

  if (!is.null(expr_matrix) && gene %in% rownames(expr_matrix)) {
    expr_gene <- expr_matrix[gene, , drop = FALSE]
    expr_df <- data.frame(
      gene = gene,
      value = as.numeric(expr_gene),
      cell = colnames(expr_matrix)
    )
    
    expr_df$datasets <- meta[expr_df$cell, "datasets"]
    expr_df$group <- cellchat_combined@idents[expr_df$cell]

    df_sub <- expr_df %>% filter(group == cell_type)

    if (length(unique(df_sub$datasets)) == 2) {
      pval <- wilcox.test(value ~ datasets, data = df_sub)$p.value
      cat("✓", gene, "-", cell_type, "- p =", signif(pval, 3), "\n")
      pval_results <- rbind(pval_results, data.frame(
        Gene = gene,
        CellType = cell_type,
        PValue = signif(pval, 5)
      ))
    } else {
      cat("⚠", gene, "-", cell_type, "- only one group present; skipped\n")
    }
  } else {
    cat("⚠", gene, "not found in expression matrix; skipped\n")
  }
}

# Save results
write.csv(pval_results, 
          file = paste0(output_path, "ligand_receptor_pvalues_only.csv"), 
          row.names = FALSE)

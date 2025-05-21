# ====================================================
# Slingshot Trajectory Analysis on Bcl-xL Seurat Object
# ====================================================

BiocManager::install("slingshot")

# Load libraries
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(scater)


# Step 1: Load the annotated Seurat object
seurat_obj <- readRDS("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl_annotated_noRBC.rds")

# Step 2: Filter out unwanted cell types
to_remove <- c("Oligodendrocytes", "Purkinje_Cells", "Interneurons", 
               "Vascular_Fibroblasts", "Macrophages_Microglia", "Endothelial")

# Use meta.data to get cell names that you want to keep
cells_to_keep <- rownames(seurat_obj@meta.data)[
  !(seurat_obj@meta.data$cell_type %in% to_remove)
]

# Subset the Seurat object using cell names
seurat_obj <- subset(seurat_obj, cells = cells_to_keep)

unique(seurat_obj@meta.data$cell_type)

# Step 3: Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_obj)

# Step 4: Run Slingshot and Save trajectory plot
sce <- slingshot(
  sce,
  clusterLabels = sce$seurat_clusters,
  reducedDim = "UMAP"
)


png(
  filename = "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/slingshot_trajectory_umap_filtered.png",
  width = 7, height = 7, units = "in", res = 300
)

plot(
  reducedDims(sce)$UMAP,
  col = sce$seurat_clusters,
  pch = 16,
  cex = 0.3,
  main = "Slingshot Trajectory (Filtered Cell Types)"
)
lines(SlingshotDataSet(sce), lwd = 2)

dev.off()


# ---------------------------------------------
# Step 5: Extract and save pseudotime values
# ---------------------------------------------
pseudotime <- slingPseudotime(sce)
head(pseudotime)

# Save pseudotime to Seurat object (first lineage)
seurat_obj$pseudotime <- pseudotime[, 1]

# Save the updated Seurat object with pseudotime added
saveRDS(
  seurat_obj,
  file = "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/bclxl_with_pseudotime_filtered.rds"
)



## Plot Pseudotime for individual Lineage
# Extract UMAP coordinates and pseudotime
umap_df <- as.data.frame(reducedDims(sce)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$pseudotime <- slingPseudotime(sce)[, 1]  # First lineage

# Plot with pseudotime gradient
library(ggplot2)

png("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/umap_pseudotime_gradient_first_lineage_filtered.png",
    width = 6, height = 5, units = "in", res = 300)

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(option = "plasma", na.value = "grey90") +
  labs(title = "UMAP Colored by Slingshot Pseudotime First lineage") +
  theme_classic()

dev.off()


## Plot Pseudotime for combined Lineage
# ------------------------------------------------
# Extract UMAP coordinates and combined pseudotime
# ------------------------------------------------
umap_df <- as.data.frame(reducedDims(sce)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

# Combine pseudotime across all lineages (minimum pseudotime per cell)
pseudotime_all <- slingPseudotime(sce)
umap_df$pseudotime <- apply(pseudotime_all, 1, function(x) min(x, na.rm = TRUE))
umap_df$pseudotime[!is.finite(umap_df$pseudotime)] <- NA  # Clean up

# ------------------------------------------------
# Plot UMAP with pseudotime gradient
# ------------------------------------------------
library(ggplot2)

png("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/umap_combined_pseudotime_gradient_filtered.png",
    width = 6, height = 5, units = "in", res = 300)

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(option = "plasma", na.value = "grey90") +
  labs(title = "UMAP Colored by Combined Pseudotime (All Lineages)") +
  theme_classic()

dev.off()


## Code to Plot individual Genes Across Pseudotime
## Loop to Plot Pseudotime for Each Lineage
# Load required libraries
library(ggplot2)
library(viridis)

# Get UMAP coordinates
umap_df_base <- as.data.frame(reducedDims(sce)$UMAP)
colnames(umap_df_base) <- c("UMAP_1", "UMAP_2")

# Extract all pseudotime values
pseudotime_all <- slingPseudotime(sce)

# Loop through each lineage
for (i in 1:ncol(pseudotime_all)) {
  
  # Clone base UMAP data and add pseudotime for this lineage
  umap_df <- umap_df_base
  umap_df$pseudotime <- pseudotime_all[, i]
  
  # Clean up NAs (those not assigned to this lineage)
  umap_df$pseudotime[!is.finite(umap_df$pseudotime)] <- NA
  
  # Set output filename
  filename <- paste0("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/umap_pseudotime_gradient_lineage_", i, "_filtered.png")
  
  # Save plot
  png(filename, width = 6, height = 5, units = "in", res = 300)
  
  print(
    ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
      geom_point(size = 0.5) +
      scale_color_viridis_c(option = "plasma", na.value = "grey90") +
      labs(title = paste("UMAP Colored by Slingshot Pseudotime - Lineage", i)) +
      theme_classic()
  )
  
  dev.off()
}



## Loop Code to Plot Multiple Genes Across Pseudotime

# Full list of genes to include (from previous and new)
genes_to_plot <- c("Mki67", "App", "Nrcam", "Gja1", "Ncam1", "Pcdhga9", "Psap", "Cadm1", "Bsg")

# Loop through each gene and generate a pseudotime plot
for (gene in genes_to_plot) {
  
  # Add expression to UMAP dataframe
  umap_df$gene_expr <- seurat_obj[["RNA"]]@data[gene, ]
  
  # Define output file path
  output_file <- paste0("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/gene_", gene, "_vs_pseudotime_filtered.png")
  
  # Save the plot
  png(output_file, width = 6, height = 4, units = "in", res = 300)
  
  # Plot expression vs pseudotime
  print(
    ggplot(umap_df, aes(x = pseudotime, y = gene_expr)) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_smooth(method = "loess", color = "firebrick") +
      labs(title = paste(gene, "Expression over Pseudotime"),
           x = "Pseudotime", y = "Log-normalized Expression") +
      theme_minimal()
  )
  
  dev.off()
}



# ============================================================
# Identify Branching Points and Diverging Fates
# ============================================================

# ----- Step 1: Check number of inferred lineages -----
n_lineages <- length(slingCurves(sce))
print(paste("Number of inferred lineages:", n_lineages))

# ----- Step 2: Plot UMAP with Slingshot trajectories -----
png("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/umap_with_trajectories_filtered.png",
    width = 6, height = 5, units = "in", res = 300)

plot(reducedDims(sce)$UMAP, col = as.factor(seurat_obj$cell_type), pch = 16, cex = 0.3,
     main = "UMAP with Slingshot Trajectories")
lines(SlingshotDataSet(sce), lwd = 2)

dev.off()

# ----- Step 3: Identify branching cells from lineage weights -----
sling_lineage_weights <- slingCurveWeights(sce)

# Cells with significant weight in more than one lineage = likely branching points
branching_cells <- which(rowSums(sling_lineage_weights > 0.25) > 1)

# Mark cells as branching or not
umap_df$branching <- ifelse(rownames(umap_df) %in% names(branching_cells), "Yes", "No")

# ----- Step 4: Plot UMAP highlighting branching points -----
png("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/umap_branching_cells_filtered.png",
    width = 6, height = 5, units = "in", res = 300)

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = branching)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("Yes" = "purple", "No" = "lightgray")) +
  labs(title = "Branching Cells in UMAP") +
  theme_classic() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
      ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

dev.off()

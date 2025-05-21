# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Define the file path
bclxl_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl_no_bclxl_all_noB6.rds"

# Load the Seurat object
bclxl <- readRDS(bclxl_path)

# Update the Seurat object
bclxl <- UpdateSeuratObject(bclxl)

# Check if the update was successful
bclxl

# View metadata structure
head(bclxl@meta.data)

table(bclxl$seurat_clusters)

# Define the file path
save_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl/UMAP_BclXL_KO_vs_Control.png"

# Create the UMAP plot
umap_plot <- DimPlot(bclxl, reduction = "umap", label = TRUE, group.by = "geno") + 
             ggtitle("Bcl-xL Knockout vs Control")


# Save the plot as a PNG file
ggsave(filename = save_path, plot = umap_plot, width = 8, height = 6, dpi = 300)


# UMAP plot colored by clusters
DimPlot(bclxl, reduction = "umap", label = TRUE, group.by = "SCT_snn_res.0.8") +
  ggtitle("Seurat Clusters (Resolution 0.8) - Bcl-xL") +
  theme_minimal()

  # Define file path to save the plot
save_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl/UMAP_BclXL_Seurat_Clusters_0.8.png"

# Save the plot
ggsave(filename = save_path, plot = last_plot(), width = 8, height = 6, dpi = 300)





# Load necessary libraries
library(Seurat)
library(SingleR)
library(celldex)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)

# Define the file path to the Seurat object
bclxl_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl_no_bclxl_all_noB6.rds"

# Load the Seurat object
bclxl <- readRDS(bclxl_path)

# Update the Seurat object
bclxl <- UpdateSeuratObject(bclxl)

# Check if the update was successful
print(bclxl)

# View metadata structure
head(bclxl@meta.data)

# Check cluster distribution
table(bclxl$seurat_clusters)

# ================================
# UMAP Visualization Before SingleR
# ================================

# Define file path to save UMAP plot (Genotype)
save_path_geno <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl/UMAP_BclXL_KO_vs_Control.png"

# Create and save the UMAP plot grouped by genotype
umap_plot_geno <- DimPlot(bclxl, reduction = "umap", label = TRUE, group.by = "geno") + 
                  ggtitle("Bcl-xL Knockout vs Control")

ggsave(filename = save_path_geno, plot = umap_plot_geno, width = 8, height = 6, dpi = 300)

# Define file path to save UMAP plot (Clusters)
save_path_clusters <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl/UMAP_BclXL_Seurat_Clusters_0.8.png"

# Create and save the UMAP plot grouped by clusters
umap_plot_clusters <- DimPlot(bclxl, reduction = "umap", label = TRUE, group.by = "SCT_snn_res.0.8") +
                      ggtitle("Seurat Clusters (Resolution 0.8) - Bcl-xL") +
                      theme_minimal()

ggsave(filename = save_path_clusters, plot = umap_plot_clusters, width = 8, height = 6, dpi = 300)

# ================================
# Perform Cell Type Annotation with SingleR
# ================================

# Extract normalized RNA counts for SingleR (NOT SCT!)
data <- GetAssayData(bclxl, assay = "RNA", slot = "data")
if (is.null(rownames(data))) {
  rownames(data) <- rownames(bclxl[["RNA"]])
}

# Load the Mouse reference dataset
ref <- celldex::MouseRNAseqData()  # Alternative: celldex::ImmGenData() for immune-focused datasets

# Convert Seurat object to SingleCellExperiment format
bclxl_sce <- as.SingleCellExperiment(bclxl)

# Perform SingleR annotation
annotations <- SingleR(test = bclxl_sce, ref = ref, labels = ref$label.main)

# Add predicted labels to the Seurat object's metadata under 'cell_type'
bclxl$cell_type <- annotations$labels[match(colnames(bclxl), rownames(annotations))]

# Check the added annotations
head(bclxl$cell_type)
table(bclxl$cell_type)

# ================================
# UMAP Visualization After SingleR
# ================================

# Define file path for cell type UMAP
save_path_celltype <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl/UMAP_BclXL_Cell_Types.png"

# Create and save the UMAP plot colored by cell types
umap_plot_celltype <- DimPlot(bclxl, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 3) +
                      ggtitle("Cell Type Annotation (SingleR-MouseRNAseqData) - Bcl-xL")


print(umap_plot_celltype)

ggsave(filename = save_path_celltype, plot = umap_plot_celltype, width = 10, height = 8, dpi = 300)

# ================================
# Save the Updated Seurat Object
# ================================
saveRDS(bclxl, "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl/bclxl_annotated.rds")



# Compare Seurat clusters and SingleR labels directly
table(bclxl$seurat_clusters, bclxl$cell_type)

# To view the table in a more interactive way in RStudio, you can use View():
View(table(bclxl$seurat_clusters, bclxl$cell_type))




# Load necessary libraries
library(Seurat)
library(SingleR)
library(celldex)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)

# Define the file path to the Seurat object
bclxl_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl_no_bclxl_all_noB6.rds"

# Load the Seurat object
bclxl <- readRDS(bclxl_path)

# Update the Seurat object
bclxl <- UpdateSeuratObject(bclxl)

# Check if the update was successful
print(bclxl)

# View metadata structure
head(bclxl@meta.data)

# Check cluster distribution
table(bclxl$seurat_clusters)

# UMAP Visualization Before SingleR

# Define file path to save UMAP plot (Genotype)
save_path_geno <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl/UMAP_BclXL_KO_vs_Control.png"

# Create and save the UMAP plot grouped by genotype
umap_plot_geno <- DimPlot(bclxl, reduction = "umap", label = TRUE, group.by = "geno") + 
                  ggtitle("Bcl-xL Knockout vs Control")

ggsave(filename = save_path_geno, plot = umap_plot_geno, width = 8, height = 6, dpi = 300)

# Define file path to save UMAP plot (Clusters)
save_path_clusters <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl/UMAP_BclXL_Seurat_Clusters_0.8.png"

# Create and save the UMAP plot grouped by clusters
umap_plot_clusters <- DimPlot(bclxl, reduction = "umap", label = TRUE, group.by = "SCT_snn_res.0.8") +
                      ggtitle("Seurat Clusters (Resolution 0.8) - Bcl-xL") +
                      theme_minimal()

ggsave(filename = save_path_clusters, plot = umap_plot_clusters, width = 8, height = 6, dpi = 300)

# Perform Cell Type Annotation with SingleR using ImmGenData()

# Extract normalized RNA counts for SingleR
data <- GetAssayData(bclxl, assay = "RNA", slot = "data")
if (is.null(rownames(data))) {
  rownames(data) <- rownames(bclxl[["RNA"]])
}

# Load the ImmGen reference dataset (better for immune-related data)
ref <- celldex::ImmGenData()

# Convert Seurat object to SingleCellExperiment format
bclxl_sce <- as.SingleCellExperiment(bclxl)

# Perform SingleR annotation
annotations <- SingleR(test = bclxl_sce, ref = ref, labels = ref$label.main)

# Add predicted labels to the Seurat object's metadata under 'cell_type'
bclxl$cell_type <- annotations$labels[match(colnames(bclxl), rownames(annotations))]

# Check the added annotations
head(bclxl$cell_type)
table(bclxl$cell_type)

# Compare Seurat clusters and SingleR labels
cluster_celltype_table <- table(bclxl$seurat_clusters, bclxl$cell_type)
View(cluster_celltype_table)

# UMAP Visualization After SingleR

# Define file path for cell type UMAP
save_path_celltype <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl/UMAP_BclXL_Cell_Types_02.png"

# Create and save the UMAP plot colored by cell types
umap_plot_celltype <- DimPlot(bclxl, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 3) +
                      ggtitle("Cell Type Annotation (SingleR - ImmGen) - Bcl-xL")

ggsave(filename = save_path_celltype, plot = umap_plot_celltype, width = 8, height = 6, dpi = 300)

# Save the Updated Seurat Object
saveRDS(bclxl, "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl/bclxl_annotated_ImmGen.rds")




# ===========================================
# Seurat Differential Gene Expression Analysis for (bclxl_no_bclxl_all_noB6.rds)
# ===========================================
# This script performs differential gene expression analysis using a Seurat object.
# It identifies differentially expressed genes, generates violin and dot plots,
# and saves the processed results.

# ===========================================
# Load Required Libraries
# ===========================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# ===========================================
# Define Paths and Helper Functions
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
fig_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl/"

# Function to Save Figures
SaveFigure <- function(plots, name, width, height, res = 200) {
  file_path <- paste0(fig_path, name, ".png")
  png(file_path, width = width, height = height, units = "in", res = res)
  print(plots)
  dev.off()
}

# ===========================================
# Load the Seurat Object
# ===========================================
seurat_filename <- "bclxl_no_bclxl_all_noB6"
bclxl <- readRDS(paste0(data_path, seurat_filename, ".rds"))

# Update the Seurat object
bclxl <- UpdateSeuratObject(bclxl)

# Check if the update was successful
print(bclxl)

# View metadata structure
head(bclxl@meta.data)

# Check cluster distribution
table(bclxl$seurat_clusters)


# ===========================================
# Normalization and Identification of Variable Features
# ===========================================
# Normalize the data
bclxl <- NormalizeData(bclxl, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
bclxl <- FindVariableFeatures(bclxl, selection.method = "vst", nfeatures = 2000)

# Identify the top 10 most highly variable genes
top10 <- head(VariableFeatures(bclxl), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(bclxl)
plot2 <- LabelPoints(plot = plot1, points = top10, xnudge = 0, ynudge = 0, repel = TRUE)
SaveFigure((plot1 + plot2), paste0(seurat_filename, "_var_features"), width = 18, height = 9)

# ===========================================
# UMAP Visualization for Clusters
# ===========================================
# Define file path for UMAP plot (clusters)
save_path_umap <- paste0(fig_path, seurat_filename, "_UMAP_clusters.png")

# Create and save the UMAP plot grouped by clusters
umap_plot <- DimPlot(bclxl, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
             ggtitle("UMAP of Seurat Clusters - Bcl-xL")

SaveFigure(umap_plot, paste0(seurat_filename, "_UMAP_clusters"), width = 10, height = 8)


# ===========================================
# Differential Gene Expression Analysis
# ===========================================
# Identify differentially expressed genes in each cluster
bclxl_markers <- FindAllMarkers(bclxl, only.pos = TRUE)

# View top markers for each cluster
bclxl_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# Save DEGs to file
saveRDS(bclxl_markers, file = paste0(data_path, seurat_filename, "_markers.RDS"))
write.csv(bclxl_markers, file = paste0(data_path, seurat_filename, "_cluster_markers.csv"), 
          sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# ===========================================
# Visualization of Key Marker Genes
# ===========================================
# Selected marker genes for visualization (unchanged from original)
selected_features <- c("Aqp4", "Pecam1", "Sox10", "Mrc1", "Dcn", "Meg3", 
                       "Ptprc", "Mki67", "Ccnd1", "Nhlh1", "Gabra2")

# ---------- Violin Plot for Key Genes ----------
plot <- VlnPlot(bclxl, features = selected_features, group.by = "seurat_clusters", pt.size = 0, ncol = 3)
SaveFigure(plot, paste0(seurat_filename, "_violin_exp1"), width = 20, height = 8)

# Violin plot with raw counts (log scale)
plot <- VlnPlot(bclxl, features = selected_features, layer = "counts", log = TRUE, group.by = "seurat_clusters", pt.size = 0, ncol = 3)
SaveFigure(plot, paste0(seurat_filename, "_violin_exp2"), width = 20, height = 8)

# ---------- Feature Plots for Selected Genes ----------
feature_plot <- FeaturePlot(bclxl, features = selected_features, ncol = 4)
SaveFigure(feature_plot, paste0(seurat_filename, "_feature_exp"), width = 20, height = 12)

# Count the number of cells in each cluster
table(bclxl$seurat_clusters)

# ---------- DotPlot for Top 5 Genes per Cluster ----------
# Identify top 5 differentially expressed genes per cluster
top5 <- bclxl_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
to_plot <- unique(top5$gene)

# Generate and save dot plot
plot <- DotPlot(bclxl, features = to_plot, group.by = "seurat_clusters") + coord_flip()
SaveFigure(plot, paste0(seurat_filename, "_dotplot_top5"), width = 9, height = 20)




# ===========================================
# Seurat Differential Gene Expression Analysis for (bclxl_clustered_res0.5.rds)
# ===========================================
# This script performs differential gene expression analysis using a Seurat object.
# It identifies differentially expressed genes, generates violin and dot plots,
# and saves the processed results.

# ===========================================
# Load Required Libraries
# ===========================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# ===========================================
# Define Paths and Helper Functions
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
fig_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl/"

# Function to Save Figures
SaveFigure <- function(plots, name, width, height, res = 200) {
  file_path <- paste0(fig_path, name, ".png")
  png(file_path, width = width, height = height, units = "in", res = res)
  print(plots)
  dev.off()
}

# ===========================================
# Load the Seurat Object
# ===========================================
seurat_filename <- "bclxl_no_bclxl_all_noB6"
bclxl <- readRDS(paste0(data_path, seurat_filename, ".rds"))

# Update the Seurat object
bclxl <- UpdateSeuratObject(bclxl)

# Check if the update was successful
print(bclxl)

# View metadata structure
head(bclxl@meta.data)

# Check cluster distribution
table(bclxl$seurat_clusters)

# Ensure the correct assay is active before clustering
DefaultAssay(bclxl) <- "SCT"

# Compute the Shared Nearest Neighbor (SNN) Graph again (only needed if not already done)
#bclxl <- FindNeighbors(bclxl, dims = 1:16)  # Use same number of dimensions as before

# Re-run clustering with the desired resolution (0.5)
bclxl <- FindClusters(bclxl, resolution = 0.5)

# View the new cluster distribution
table(bclxl$SCT_snn_res.0.5)

# Generate a UMAP plot colored by the new clusters
umap_plot <- DimPlot(bclxl, reduction = "umap", label = TRUE, group.by = "SCT_snn_res.0.5") +
             ggtitle("UMAP of Seurat Clusters - Resolution 0.5")

# Save the plot
SaveFigure(umap_plot, "UMAP_Seurat_Clusters_0.5", width = 10, height = 8)

# Save the updated Seurat object
saveRDS(bclxl, file = paste0(data_path, "bclxl_clustered_res0.5.rds"))


# ===========================================
# Load the Updated Seurat Object
# ===========================================
seurat_filename <- "bclxl_clustered_res0.5"
bclxl <- readRDS(paste0(data_path, seurat_filename, ".rds"))

# Set SCT as the active assay for feature selection and analysis
DefaultAssay(bclxl) <- "SCT"

# ===========================================
# Identification of Highly Variable Features
# ===========================================
# Identify highly variable features (no need to normalize again)
bclxl <- FindVariableFeatures(bclxl, selection.method = "vst", nfeatures = 2000)

# Identify the top 10 most highly variable genes
top10 <- head(VariableFeatures(bclxl), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(bclxl)
plot2 <- LabelPoints(plot = plot1, points = top10, xnudge = 0, ynudge = 0, repel = TRUE)
SaveFigure((plot1 + plot2), paste0(seurat_filename, "_var_features"), width = 18, height = 9)

# ===========================================
# UMAP Visualization for Clusters
# ===========================================
# Create and save the UMAP plot grouped by clusters at resolution 0.5
umap_plot <- DimPlot(bclxl, reduction = "umap", label = TRUE, group.by = "SCT_snn_res.0.5") +
             ggtitle("UMAP of Seurat Clusters - Resolution 0.5")

SaveFigure(umap_plot, paste0(seurat_filename, "_UMAP_clusters_res0.5"), width = 10, height = 8)

# ===========================================
# Differential Gene Expression Analysis
# ===========================================
# Identify differentially expressed genes in each cluster
bclxl_markers <- FindAllMarkers(bclxl, only.pos = TRUE)

# View top markers for each cluster
bclxl_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# Save DEGs to file
saveRDS(bclxl_markers, file = paste0(data_path, seurat_filename, "_markers.RDS"))
write.csv(bclxl_markers, file = paste0(data_path, seurat_filename, "_cluster_markers.csv"), 
          sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# ===========================================
# Visualization of Key Marker Genes
# ===========================================
# Selected marker genes for visualization (unchanged from original)
selected_features <- c("Aqp4", "Pecam1", "Sox10", "Mrc1", "Dcn", "Meg3", 
                       "Ptprc", "Mki67", "Ccnd1", "Nhlh1", "Gabra2")

# ---------- Violin Plot for Key Genes ----------
plot <- VlnPlot(bclxl, features = selected_features, group.by = "SCT_snn_res.0.5", pt.size = 0, ncol = 3)
SaveFigure(plot, paste0(seurat_filename, "_violin_exp1"), width = 20, height = 8)

# Violin plot with raw counts (log scale)
plot <- VlnPlot(bclxl, features = selected_features, layer = "counts", log = TRUE, group.by = "SCT_snn_res.0.5", pt.size = 0, ncol = 3)
SaveFigure(plot, paste0(seurat_filename, "_violin_exp2"), width = 20, height = 8)

# ---------- Feature Plots for Selected Genes ----------
feature_plot <- FeaturePlot(bclxl, features = selected_features, ncol = 4)
SaveFigure(feature_plot, paste0(seurat_filename, "_feature_exp"), width = 20, height = 12)

# Count the number of cells in each cluster
table(bclxl$SCT_snn_res.0.5)

# ---------- DotPlot for Top 5 Genes per Cluster ----------
# Identify top 5 differentially expressed genes per cluster
top5 <- bclxl_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
to_plot <- unique(top5$gene)

# Generate and save dot plot
plot <- DotPlot(bclxl, features = to_plot, group.by = "SCT_snn_res.0.5") + coord_flip()
SaveFigure(plot, paste0(seurat_filename, "_dotplot_top5"), width = 9, height = 20)




## annotate cell type

cluster_annotations <- c(
  "0"  = "Stem_A",
  "1"  = "Mitotic_A",
  "2"  = "Early_Differentiating",
  "3"  = "Stem_B",
  "4"  = "Late_Differentiating",
  "5"  = "Mitotic_B",
  "6"  = "Oligodendrocytes",
  "7"  = "Red_Blood_Cells",
  "8"  = "Purkinje_Cells",
  "9"  = "Macrophages_Microglia",
  "10" = "Astrocytes",
  "11" = "Pericytes",
  "12" = "Interneurons",
  "13" = "Vascular_Fibroblasts",
  "14" = "Endothelial",
  "15" = "Differentiated_CGN_like"
)

# Confirm Unique Cluster Labels
unique(bclxl$SCT_snn_res.0.5)

# If the output contains levels, convert to character
bclxl$SCT_snn_res.0.5 <- as.character(bclxl$SCT_snn_res.0.5)

# Check again to ensure conversion worked
unique(bclxl$SCT_snn_res.0.5)

# Ensure cluster_annotations names are explicitly characters
names(cluster_annotations) <- as.character(names(cluster_annotations))

# Assign Cell Type Labels to the Seurat Object
bclxl$cell_type <- unname(cluster_annotations[bclxl$SCT_snn_res.0.5])

# Verify Cluster-to-Cell Type Mapping
table(bclxl$cell_type)

# Verify the Annotation
table(bclxl$SCT_snn_res.0.5, bclxl$cell_type)

# Check if Any NA Values Exist
sum(is.na(bclxl$cell_type))

# UMAP plot with annotated cell types
umap_annotated <- DimPlot(bclxl, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 4) +
                  ggtitle("UMAP of Annotated Cell Types - Bcl-xL")

# Save the annotated UMAP plot
SaveFigure(umap_annotated, paste0(seurat_filename, "_UMAP_celltypes"), width = 10, height = 8)

saveRDS(bclxl, file = paste0(data_path, "bclxl_annotated.rds"))




## Remove cluster "7" (Red Blood Cells) from the annotations
# ===========================================
# Define Paths and Helper Functions
# ===========================================
data_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/"
fig_path <- "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl/"

# Define the base filename for saving figures
seurat_filename <- "bclxl_annotated_noRBC"

# Function to Save Figures
SaveFigure <- function(plots, name, width, height, res = 200) {
  file_path <- paste0(fig_path, name, ".png")
  png(file_path, width = width, height = height, units = "in", res = res)
  print(plots)
  dev.off()
}

# Load the previously saved Seurat object
bclxl <- readRDS(file = paste0(data_path, "bclxl_annotated.rds"))

# Remove cells belonging to cluster "7" (Red Blood Cells)
bclxl <- subset(bclxl, subset = SCT_snn_res.0.5 != "7")

# Define updated cluster annotations (without Red Blood Cells)
cluster_annotations <- c(
  "0"  = "Stem_A",
  "1"  = "Mitotic_A",
  "2"  = "Early_Differentiating",
  "3"  = "Stem_B",
  "4"  = "Late_Differentiating",
  "5"  = "Mitotic_B",
  "6"  = "Oligodendrocytes",
  "8"  = "Purkinje_Cells",
  "9"  = "Macrophages_Microglia",
  "10" = "Astrocytes",
  "11" = "Pericytes",
  "12" = "Interneurons",
  "13" = "Vascular_Fibroblasts",
  "14" = "Endothelial",
  "15" = "Differentiated_CGN_like"
)

# Ensure cluster_annotations names are explicitly characters
names(cluster_annotations) <- as.character(names(cluster_annotations))

# Assign Cell Type Labels to the Seurat Object
bclxl$cell_type <- unname(cluster_annotations[bclxl$SCT_snn_res.0.5])

# Verify Cluster-to-Cell Type Mapping
table(bclxl$cell_type)

# UMAP plot with updated annotated cell types
umap_annotated <- DimPlot(bclxl, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 4) +
                  ggtitle("Updated UMAP of Annotated Cell Types - Bcl-xL (No RBC)")

# Save the updated UMAP plot
SaveFigure(umap_annotated, paste0(seurat_filename, "_UMAP_celltypes_updated"), width = 10, height = 8)

# Save the updated Seurat object (without cluster 7)
saveRDS(bclxl, file = paste0(data_path, "bclxl_annotated_noRBC.rds"))


# make sure Red_Blood_Cells have been removed
# Load the previously saved Seurat object
bclxl <- readRDS(file = paste0(data_path, "bclxl_annotated_noRBC.rds"))

bclxl <- readRDS(file = paste0(data_path, "bclxl_annotated.rds"))

head(bclxl@meta.data)

table(bclxl@meta.data$cell_type)

bclxl


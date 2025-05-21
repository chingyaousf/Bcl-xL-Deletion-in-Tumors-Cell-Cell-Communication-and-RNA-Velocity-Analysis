# Bcl-xL Deletion in Tumors: Cell-Cell Communication and RNA Velocity Analysis

## Background:

Bcl-xL is an anti-apoptotic protein known for its role in promoting tumor survival. To investigate how Bcl-xL deletion alters tumor cell dynamics and intercellular signaling, we conducted a comprehensive single-cell RNA-sequencing analysis comparing Bcl-xL knockout (KO) and heterozygous (Het) tumor samples. This project integrates cell type annotation, intercellular communication inference using **CellChat**, and lineage trajectory analysis using **Slingshot**.

**CellChat** is a computational framework designed to infer and analyze intercellular communication networks from single-cell transcriptomic data. In this study, we used CellChat to identify signaling interactions that are enriched or diminished following Bcl-xL deletion. This allowed us to systematically assess how cell-cell communication is rewired in the tumor microenvironment.

**Slingshot** is a trajectory inference method that models cellular differentiation by estimating pseudotime and lineage relationships from low-dimensional embeddings. We employed Slingshot to examine how Bcl-xL deletion influences the developmental trajectories of tumor-associated cells, revealing shifts in differentiation patterns that may be linked to altered signaling environments.

## Methods:

### Data Preprocessing and Annotation

-   Seurat objects were loaded and preprocessed using standard workflows including normalization, dimensionality reduction, and clustering.

-   Cell types were annotated based on an existing reference list of previously identified and validated cell types.

-   Red blood cells were excluded prior to downstream analyses to improve interpretability.

### CellChat Analysis

-   CellChat was used to infer intercellular communication networks separately for KO and Het samples.

-   Network metrics such as interaction strength and number were computed.

-   Signaling pathways were analyzed to identify ligand-receptor activity changes.

-   Differential interaction networks were visualized using heatmaps, scatter plots, and chord diagrams.

-   Specific attention was given to pathways with strong structural differences, such as CNTN, GAP, NCAM, PCDH, and PSAP.

### Trajectory Inference

-   Pseudotime and differentiation trajectories were estimated using Slingshot on UMAP embeddings.

-   Both raw and filtered datasets (excluding non-relevant cell types) were processed.

-   Pseudotime gradients were visualized to assess progression across lineages.

-   Expression of specific genes across pseudotime was plotted to examine dynamic transcriptional changes.

## Results:

-   Fifteen distinct cell types were annotated post-filtering, including astrocytes, stem-like populations, and mitotic subtypes.

-   KO tumors showed altered communication profiles compared to Het controls, with significant rewiring in neural and adhesion-related signaling pathways.

-   Several ligand-receptor interactions were upregulated in KO, including those involved in structural cell-cell interactions.

-   Trajectory inference revealed altered developmental paths in KO tumors, indicating changes in cellular differentiation states.

-   Combined analyses from CellChat and Slingshot suggest that Bcl-xL deletion not only reshapes the intercellular communication landscape but also modulates the progression and identity of tumor cell populations.

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/bclxl/bclxl_clustered_res0.5_UMAP_celltypes.png?raw=true){width="500"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/bclxl/bclxl_clustered_res0.5_var_features.png?raw=true){width="500"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/bclxl/bclxl_clustered_res0.5_feature_exp.png?raw=true)

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/bclxl/bclxl_clustered_res0.5_violin_exp2.png?raw=true)

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/bclxl/bclxl_clustered_res0.5_dotplot_top5.png?raw=true)

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/cellchat_het_02/interaction_number_het_02.png?raw=true){width="500"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/cellchat_het_02/interaction_strength_het_02.png?raw=true){width="500"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/cellchat_ko_02/interaction_number_KO_02.png?raw=true){width="500"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/cellchat_ko_02/interaction_strength_KO_02.png?raw=true){width="500"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/diffInteraction_het_ko.png?raw=true){width="500"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/diffStrength_het_ko.png?raw=true){width="500"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/differential_interactions_heatmap_het_ko.png?raw=true){width="1000"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/signalingRole_scatter_comparison_het_ko.png?raw=true)

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/scatter_incoming_centrality_with_count.png?raw=true){width="700"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/informationFlow_structuralPathways_het_ko.png?raw=true){width="700"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/informationFlow_structuralPathways_gg2_het_ko.png?raw=true){width="700"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/compare_interactions_het_ko.png?raw=true){width="500"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/APP_chord_comparison_het_ko.png?raw=true){width="700"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/file/violin_LR_comparison_APP_het_ko.png?raw=true){width="700"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/PCDH_chord_comparison_het_ko.png?raw=true){width="700"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/file/violin_LR_comparison_PCDH_het_ko.png?raw=true){width="700"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/CADM_chord_comparison_het_ko.png?raw=true){width="700"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/file/violin_LR_comparison_CADM_het_ko.png?raw=true){width="700"}

![](https://github.com/chingyaousf/Bcl-xL-Deletion-in-Tumors-Cell-Cell-Communication-and-RNA-Velocity-Analysis/blob/main/plots/comparasion_het_ko_02/distance_structural_similarity_het_ko_02.png?raw=true){width="800"}

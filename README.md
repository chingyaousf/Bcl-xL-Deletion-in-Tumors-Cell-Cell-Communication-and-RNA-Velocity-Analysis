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

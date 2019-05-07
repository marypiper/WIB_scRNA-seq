---
title: "Single-cell RNA-seq: Clustering Analysis"
author: "Mary Piper, Lorena Pantano, Meeta Mistry, Radhika Khetani"
date: Tuesday, April 24, 2019
---

Approximate time: 90 minutes

## Learning Objectives:

* Understand how to determine the most variable genes for the gene expression
* Utilize methods for evaluating the selection of PCs to use for clustering
* Perform clustering of cells based on significant PCs
* Evaluate whether clustering artifacts are present 
* Determine the quality of clustering with PCA and tSNE plots and understand when to re-cluster

# Single-cell RNA-seq clustering analysis


Now that we have our high quality cells, we want to know the different cell types present within our population of cells. 

<img src="../img/sc_workflow.png" width="800">

_**Goals:**_ 
 
 - _To **generate cell type-specific clusters** and use known markers to determine identity of some/all of the clusters._
 - _To **determine whether clusters represent true cell types or cluster due to biological or technical variation**, such as clusters of cells in the S phase of the cell cycle, clusters of specific batches, or cells with high mitochondrial content._

_**Challenges:**_
 
 - _Clustering so that **cells of the same cell type from different conditions cluster together**_
 - _**Removing unwanted variation** so that we do not have cells clustering by artifacts_
 - _**Identifying the cell types** of each cluster_
 - _Maintaining patience as this can be a highly interative process between clustering and marker identification (sometimes even going back to the QC filtering_

_**Recommendations:**_
 
 - _Have a good idea of your expectations for the **cell types to be present** prior to performing the clustering. Know whether you expect cell types of low complexity or higher mitochondrial content or whether the cells are differentiating_
 - _If you have **more than one condition**, it's often helpful to perform integration to align the cells_
 - _**Regress out** number of UMIs, mitochondrial content, and cell cycle, if needed and appropriate for experiment, so not to drive clustering_
 - _Identify any junk clusters for removal. Possible junk clusters could include those with high **mitochondrial content** and low UMIs/genes after regression of mitochondrial content_
 - _If **not detecting all cell types as separate clusters**, try changing the resolution or the number of PCs used for clustering_


## Clustering workflow

To do this we are going to perform a clustering analysis. The workflow for this analysis is adapted from the following sources:

- Satija Lab: [Seurat v3 Guided Integration Tutorial](https://satijalab.org/seurat/v3.0/immune_alignment.html)
- Paul Hoffman: [Cell-Cycle Scoring and Regression](http://satijalab.org/seurat/cell_cycle_vignette.html)

To identify clusters, the following steps will be performed:

1. **Normalization** and **identification of high variance genes** in each sample
2. **Integration** of the samples using shared highly variable genes (optional, but recommended to align cells from different samples)
3. **Scaling** and **regression** of sources of unwanted variation (e.g. number of UMIs per cell, mitochondrial transcript abundance, cell cycle phase)
4. **Clustering cells** based on top PCs (metagenes)
5. Exploration of **quality control metrics**: determine whether clusters unbalanced wrt UMIs, genes, cell cycle, mitochondrial content, samples, etc.
6. Searching for expected cell types using **known markers**
7. **Marker identification** for each cluster


## **Normalization** and **identification of high variance genes** in each sample

The first step in the analysis is to normalize the raw counts to account for differences in sequencing depth per cell **for each sample**. By default, the raw counts are normalized using global-scaling normalization by performing the following:

1. normalizing the gene expression measurements for each cell by the total expression 
2. multiplying this by a scale factor (10,000 by default)
3. log-transforming the result

Following normalization, we want to **identify the most variable genes** (highly expressed in some cells and lowly expressed in others) to use for downstream clustering analyses. 

The mean-variance relationship of the data is modelled, and the 2,000 most variable genes are returned.
 
<p align="center">
<img src="../img/variable_genes.png" width="400">
</p>

_Image credit: [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html)_

> **NOTE:** Seurat has just incorporated the `sctransform` tool for better normalization, scaling, and finding of variable genes. There is a new [vignette](https://satijalab.org/seurat/v3.0/sctransform_vignette.html) and [preprint](https://www.biorxiv.org/content/biorxiv/early/2019/03/18/576827.full.pdf) available to explore this new methodology.

At this point in the workflow, we can either:

- use the most variable genes to identify significant principal components to use for clustering (single condition analyses)
- use the most variable genes that are shared between samples to correlate the components for integration of the samples (multiple condition analyses or batches)

Oftentimes, if samples are created using different conditions (or batches), the effect of the condition can be so great that the cells will cluster by condition instead of cell type. To ensure the same cell types cluster together, we can perform an integration analysis detailed in the paper by [Stuart and Bulter et al. (2018)](https://www.biorxiv.org/content/early/2018/11/02/460147). 

## **Integration** of the samples using shared highly variable genes (optional)

_**This step can greatly improve your clustering when you have multiple samples**. It can help to first run samples individually if unsure what clusters to expect, but when clustering the cells from multiple conditions, integration can help ensure the same cell types cluster together._

Using the shared highly variable genes from each sample, we integrate the samples to overlay cells that are similar or have a "common set of biological features". The process of integration uses canonical correlation analysis (CCA) and mutual nearest neighbors (MNN) to identify shared subpopulations across samples or datasets [[Stuart and Bulter et al. (2018)](https://www.biorxiv.org/content/early/2018/11/02/460147)]. 

Specifically, this method expects "correspondences" or **shared biological states** among at least a subset of single cells across the samples. 

1. CCA uses **shared highly variable genes** to reduce the dimensionality of the data and align the cells in each sample into the maximally correlated space (based on sets of genes exhibiting robust correlation in expression).
2. Identify mutual nearest neighbors, or 'anchors' across datasets (sometimes incorrect anchors are identified)
3. Assess the similarity between anchor pairs by the overlap in their local neighborhoods (incorrect anchors will have low scores)
4. Use anchors and corresponding scores to transform cell expression values, allowing for the integration of the datasets (different samples, datasets, modalities)
	- Transformation of each cell uses a weighted average of the two cells of each anchor across anchors of the datasets. Weights determined by cell similarity score (distance between cell and k nearest anchors) and anchor scores, so cells in the same neighborhood should have similar correction values.

If cell types are present in one dataset, but not the other, then the cells will still appear as a separate sample-specific cluster.

<p align="center">
<img src="../img/integration.png" width="600">
</p>

_Image credit: Stuart T and Butler A, et al. Comprehensive integration of single cell data, bioRxiv 2018 (https://doi.org/10.1101/460147)_

## Scaling and regression of sources of unwanted variation

In addition to the interesting variation in your dataset that separates the different cell types, there is also "uninteresting" sources of variation present that can obscure the cell type-specific differences. This can include technical noise, batch effects, and/or uncontrolled biological variation (e.g. cell cycle).

### Cell cycle scoring

Cell cycle variation is a common source of uninteresting variation in single-cell RNA-seq data. To examine cell cycle variation in our data, we assign each cell a score, based on its expression of G2/M and S phase markers. 

> An overview of the cell cycle phases is given in the image below:
> 
> <p align="center">
><img src="../img/cell_cycle.png" width="200">
></p> 	
> 
> _Adapted from [Wikipedia](https://en.wikipedia.org/wiki/Cell_cycle) (Image License is [CC BY-SA 3.0](https://en.wikipedia.org/wiki/Wikipedia:Text_of_Creative_Commons_Attribution-ShareAlike_3.0_Unported_License))_
> 
> - **G0:** Quiescence or resting phase. The cell is not actively dividing, which is common for cells that are fully differentiated. Some types of cells enter G0 for long periods of time (many neuronal cells), while other cell types never enter G0 by continuously dividing (epithelial cells).
> - **G1:** Gap 1 phase represents the **beginning of interphase**. During G1 there is growth of the non-chromosomal components of the cells. From this phase, the cell may enter G0 or S phase.
> - **S:** Synthesis phase for the replication of the chromosomes (also part of interphase).
> - **G2:** Gap 2 phase represents the **end of interphase**, prior to entering the mitotic phase. During this phase th cell grows in preparation for mitosis and the spindle forms.
> - **M:** M phase is the nuclear division of the cell (consisting of prophase, metaphase, anaphase and telophase).
	

The [Cell-Cycle Scoring and Regression tutorial](https://satijalab.org/seurat/v3.0/cell_cycle_vignette.html) from Seurat makes available a list of cell cycle phase marker genes for humans, while the HBC core has [compiled lists](https://github.com/hbc/tinyatlas/tree/master/cell_cycle) for other organisms.

After scoring each gene for cell cycle phase, we can perform PCA using the expression of cell cycle genes. If the cells group by cell cycle in the PCA, then we would want to regress out cell cycle variation, **unless cells are differentiating**.  

<p align="center">
<img src="../img/SC_preregressed_phase_pca.png" width="400">
</p>


> **NOTE:** If cells are known to be differentiating and there is clear clustering differences between G2M and S phases, then you may want to regress out by the difference between the G2M and S phase scores as described in the [Seurat tutorial](https://satijalab.org/seurat/v3.0/cell_cycle_vignette.html), thereby still differentiating the cycling from the non-cycling cells.

### Apply regression variables

**Regressing variation due to uninteresting sources can improve downstream identification of principal components and clustering.** To mitigate the effects of these signals, Seurat constructs linear models to predict gene expression based on the variables to regress.

We generally recommend regressing out **number of UMIs, mitochondrial ratio, and possibly cell cycle** if needed, as a standard first-pass approach. However, if the differences in mitochondrial gene expression represent a biological phenomenon that may help to distinguish cell clusters, then we advise not regressing the mitochondrial expression.

When regressing out the effects of cell-cycle variation, include S-phase score and G2M-phase score for regression. Cell-cycle regression is generally recommended but should be avoided for samples containing cells undergoing differentiation.

> **NOTE:** If using the `sctransform` tool, there is no need to regress out number of UMIs as it is corrected for in the function.

## Clustering cells based on top PCs (metagenes)

### Identify significant PCs

To overcome the extensive technical noise in any single gene for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a "metagene" that combines information across a correlated gene set. Determining how many PCs to include downstream is therefore an important step. Often it is useful to explore the PCs prior to identifying the significant principal components to include for the downstream clustering analysis.

We can also explore the expression of the top most variant genes for select PCs using a heatmap with the genes and cells ordered by PC scores:

<img src="../img/SC_pc_heatmap.png" width="700">

However, the main analysis used to determine how many PCs to use for the downstream analysis is done through plotting the standard deviation of each PC as an elbow plot. 
Where the elbow appears is usually the threshold for identifying the most significant PCs to include.

<p align="center">
<img src="../img/SC_elbowplot.png" width="500">
</p>

While this gives us a good idea of the number of PCs to include, a more quantitative approach may be a bit more reliable.

PC selection — identifying the true dimensionality of a dataset — is an important step for our clustering analysis, but can be challenging/uncertain. While there are a variety of ways to choose a threshold, we're going to calculate where the principal components start to elbow by taking the larger value of:

1. The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
2. The point where the percent change in variation between the consequtive PCs is less than 0.1%.

### Cluster the cells

We can now use these significant PCs identified  to determine which cells exhibit similar expression patterns for clustering. To do this, Seurat uses a graph-based clustering approach, which embeds cells in a graph structure, using a K-nearest neighbor (KNN) graph (by default), with edges drawn between cells with similar gene expression patterns. Then, it attempts to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’ [[Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html)].

The `resolution` is an important argument that sets the "granularity" of the downstream clustering and will need to be optimized to the experiment.  For datasets of 3,000 - 5,000 cells, the `resolution` set between `0.4`-`1.4` generally yields good clustering. Increased resolution values lead to a greater number of clusters, which is often required for larger datasets. 

We often provide a series of resolution options during clustering, which can be used downstream to choose the best resolution.

## Creating non-linear dimensional reduction (UMAP/tSNE)

To visualize the clusters, there are a few different options that can be helpful, including t-distributed stochastic neighbor embedding (t-SNE), Uniform Manifold Approximation and Projection (UMAP), and PCA. The goals of these methods is to have similar cells closer together in low-dimensional space.

For t-SNE, the cells within the graph-based clusters just determined should co-localize on the t-SNE plot. This is because the t-SNE aims to place cells with similar local neighborhoods in high-dimensional space together in low-dimensional space. **Note that distance between clusters on the t-SNE plots does not represent degree of similarity between clusters.**

UMAP is a dimensionality reduction technique that is similar to t-SNE, but where the distances between cells represent similarity in expression. While more informative about distances, sometimes the UMAP plot can make the visualization of different clusters less spread out and differences between clusters can be hard to see. Therefore we often run both UMAP and tSNE.

<p align="center">
<img src="../img/SC_dimplot_tsne.png" width="600">
</p>

To explore similarity in gene expression between the clusters, plotting the clusters in PCA space can be more informative.

<p align="center">
<img src="../img/SC_dimplot_pca.png" width="600">
</p>


# Exploration of quality control metrics

To determine whether our clusters might be due to artifacts such as cell cycle phase or mitochondrial expression, it can be useful to explore these metrics visually to see if any clusters exhibit enrichment or are different from the other clusters. However, if enrichment or differences are observed for particular clusters it may not be worrisome if it can be explained by the cell type. 

First, we can explore cell cycle by cluster. We would also want to see if we have sample-specific clusters.

<p align="center">
<img src="../img/SC_phase_tsne_pca.png" width="600">
</p>

Next, we will explore additional metrics, such as the number of UMIs and genes per cell, S-phase and G2M-phase markers, and mitochondrial gene expression:

<img src="../img/SC_metrics_tsne.png" width="600">

We can also explore how well our clusters separate by the different PCs; we hope that the defined PCs separate the cell types well. In the tSNE plots below, the cells are colored by their PC score for each respective principal component.

<p align="center">
<img src="../img/SC_clusters_by_pc.png" width="600">
</p>

# Evaluating clustering

To determine whether our clustering and resolution are appropriate for our experiment, it is helpful to explore a handful of markers for each of the major cell types that we expect to be in our data and see how they segregate.

For example if we were interested in exploring known immune cell markers, such as:

|Marker| Cell Type|
|:---:|:---:|
|IL7R	|CD4 T cells|
|CD14, LYZ|	CD14+ Monocytes|
|MS4A1|	B cells|
|CD8A|	CD8 T cells|
|FCGR3A, MS4A7	|FCGR3A+ Monocytes|
|GNLY, NKG7|	NK cells|
|FCER1A, CST3	|Dendritic Cells|
|PPBP|	Megakaryocytes|


We can plot the expression of these genes with darker blue representing higher levels of expression. 


<img src="../img/SC_custom_genes_tsne.png" width="600">

Based on these markers, we can conjecture the identity of each of the clusters based on the canonical cell type markers:

| Cluster|Marker| Cell Type|
|:---:|:---:|:---:|
| 0-1|IL7R	|CD4 T cells|
| 2|CD14, LYZ|	CD14+ Monocytes|
| 3|MS4A1|	B cells|
| 4|CD8A|	CD8 T cells|
| 5|FCGR3A, MS4A7	|FCGR3A+ Monocytes|
| 6|GNLY, NKG7|	NK cells|
| Unidentified |FCER1A, CST3	|Dendritic Cells|
| Unidentified |PPBP|	Megakaryocytes|

Based on these results, it indicates that there are some clusters that we are not identifying that appear to be separate cell types. The megakaryocytes and the dendritic cells appear clustered with other cell type clusters, so what do we do with them? 

We would generally want to go back through the clustering, but change parameters. Did we use too few principal components that we are just not separating out these cells? We can look at our PC gene expression overlapping the tSNE plots and see these cell populations separate by PC6 and PC8, so the variation seems to be captured by our PCs. However, we might not have had a high enough resolution for our tSNE when we performed the clustering. We would want to try to re-run the tSNE with higher resolution. 

After we have identified our desired clusters, we can move on to marker identification, which will allow us to verify the identity of certain clusters and help surmise the identity of any unknown clusters. Since we have two clusters identified as CD4 T cells, we may also want to know which genes are differentially expressed between these two clusters.

[Click here for next lesson](https://github.com/marypiper/WIB_scRNA-seq/blob/master/lessons/SC_marker_identification.md)

***


*This lesson has been developed by Mary Piper. These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson were derived from work that is Copyright © [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). 
All HBC instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
* *Adapted from the lessons by [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/) with contributions from Mary Piper, Radhika Khetani, Meeta Mistry, Rory Kirchner, and Lorena Pantano.*
* *A portion of these materials and hands-on activities were adapted from the [Satija Lab's](https://satijalab.org/) [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html)*

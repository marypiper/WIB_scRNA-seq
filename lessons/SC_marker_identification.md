---
title: "Single-cell RNA-seq: Marker identification"
author: "Mary Piper, Lorena Pantano, Meeta Mistry, Radhika Khetani"
date: Tuesday, April 24, 2019
---

Approximate time: 45 minutes

## Learning Objectives:

* Understand how to determine markers of individual clusters
* Understand the iterative processes of clustering and marker identification

# Single-cell RNA-seq marker identification

Now that we have the single cells clustered based on different cell types,  we are ready to move forward with identifying cluster markers. 

<img src="../img/sc_workflow.png" width="800">

_**Goals:**_ 
 
 - _To **determine the gene markers** for each of the clusters_
 - _To **identify cell types** of each cluster using markers_
 - _To determine whether need to **re-cluster based on cell type markers**, perhaps clusters need to be merged or split

_**Challenges:**_
 
 - _Over-interpretation of the results_
 - _Combining different types of marker identification_

_**Recommendations:**_
 
 - _Identify all markers conserved between conditions for each cluster_
 - _Identify markers that are differentially expressed between specific clusters_


There are a few different types of marker identification that we can explore using Seurat. Each with their own benefits and drawbacks:

1. **Identification of all markers for each cluster:** this analysis compares each cluster against all others and outputs the genes that are differentially expressed/present. 
2. **Identification of conserved markers for each cluster regardless of condition:** This analysis looks for those genes that are conserved in the cluster across all conditions. This analysis will output genes that are consistently differentially expressed/present for all of the sample groups. These genes can help to figure out the identity for the cluster. Often, this analysis is performed only for those clusters whose identity is uncertain or novel.
3. **Marker identification between specific clusters:** this analysis explores differentially expressed genes between specific clusters. This analysis is most useful for determining differences in gene expression between clusters with markers that are similar in the above analyses. 

## Identification of all markers for each cluster

For this analysis we are comparing each cluster against all other clusters to identify cluster markers. 

To be identified as a marker, we specified that a gene needed to be detected at a minimum percentage of 0.25 in either of the two groups of cells and/or difference in expression is at least 0.25 between the two groups and/or log2 fold change is greater than 0.25.

**Usually the top markers are relatively trustworthy, but because of inflated p-values, many of the less significant genes are not so trustworthy as markers.**

When looking at the output, we suggest looking for markers with large differences in expression between `pct.1` and `pct.2` and larger fold changes. For instance if `pct.1` = 0.90 and `pct.2` = 0.80, I might not be as excited about that marker. However, if `pct.2` = 0.1 instead, then I would be much more excited about it. Also, I look for the majority of cells expressing marker in my cluster of interest. If `pct.1` is low, such as 0.3, I again might not be as interested in it.

- **cluster:** number corresponding to cluster
- **gene:** gene id
- **avg_logFC:** average log2 fold change. Positive values indicate that the gene is more highly expressed in the cluster.
- **pct.1**: The percentage of cells where the gene is detected in the cluster
- **pct.2**: The percentage of cells where the gene is detected on average in the other clusters
- **p_val:** p-value not adjusted for multiple test correction
- **p_val_adj:** Adjusted p-value, based on bonferroni correction using all genes in the dataset, used to determine significance

<img src="../img/all_markers1.png" width="750">

## Interpretation of the marker results

Using Seurat for marker identification is a rather quick and dirty way to identify markers. Usually the top markers are relatively trustworthy; however, because of inflated p-values, many of the less significant genes are not so trustworthy as markers. 

When looking at the output, we suggest looking for marker genes with large differences in expression between `pct.1` and `pct.2` and larger fold changes. For instance if `pct.1` = 0.90 and `pct.2` = 0.80 and had lower log2 fold changes, that marker might not be as exciting. However, if `pct.2` = 0.1 instead, then it would be a lot more exciting. 

<img src="../img/all_markers3.png" width="750">

If there were any clusters whose identity we were unsure of, we could looked for conserved markers.

# Assigning cell type identity to clusters

We can often go through the top markers to identify the cell types. We have to use what we know about the biology of the expected cells to determine the cell populations represented by each cluster. 

<img src="../img/tSNE.png" width="600">

To get a better idea of cell type identity we can explore the expression of different identified markers by cluster. For example, we can look at the cluster 3 markers by cluster by tSNE or violoin plot:

<img src="../img/tSNE-multiple.png" width="600">
       

<img src="../img/violinplot.png" width="600">

These results and plots can help us determine the identity of these clusters or verify what we hypothesize the identity to be after exploring the canonical markers of expected cell types previously.

**Marker identification between specific clusters:**

Sometimes the list of markers returned don't sufficiently separate some of the clusters. For instance, we had previously identified clusters 0 and 1 as CD4+ T cells, and if we would like to determine the genes that are differentially expressed between these specific clusters. 

<img src="../img/t-cell_markers2.png" width="750">

When looking through the results, the most significant marker is `ENSG00000196154`, which corresponds to **S100A4**, a gene exclusively expressed by memory T cells of CD4+ or CD8+ subpopulations. Other markers listed also indicate that cluster 0 represents naive CD4+ T cells, while cluster 1 represents memory CD4+ T cells.

Now taking all of this information, we can surmise the cell types of the different clusters. 

<img src="../img/tSNE-labelled3.png" width="600">

***


*This lesson has been developed by Mary Piper. These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson were derived from work that is Copyright © [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). 
All HBC instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
* *Adapted from the lessons by [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/) with contributions from Mary Piper, Radhika Khetani, Meeta Mistry, Rory Kirchner, and Lorena Pantano.*
* *A portion of these materials and hands-on activities were adapted from the [Satija Lab's](https://satijalab.org/) [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html)*

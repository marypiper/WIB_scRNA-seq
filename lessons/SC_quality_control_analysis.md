---
title: "Single-cell RNA-seq: Quality Control Analysis"
author: "Mary Piper, Lorena Pantano, Meeta Mistry, Radhika Khetani"
date: Monday, April 23, 2019
---

Approximate time: 90 minutes

## Learning Objectives:

* Understand how to bring in data from single-cell RNA-seq experiments
* Construct QC metrics to explore visually
* Evaluate the QC metrics to filter out low quality cells

# Single-cell RNA-seq: Quality control

The next step in the analysis is the quality control of the data at both the sample and cellular levels. This is a critical step in our analysis workflow, since failure to remove poor quality cells or a poor quality sample could make interpretation of our data more difficult.
 
<p align="center">
<img src="../img/sc_workflow.png" width="800">
</p>

_**Goals:**_ 
 
 - _To **filter the data to only include true cells that are of high quality**, so that when we cluster our cells it is easier to identify distinct cell type populations_
 - _To **identify any failed samples** and either try to salvage the data or remove from analysis, in addition to, trying to understand why the sample failed_

_**Challenges:**_
 
 - _Delineating cells that are **poor quality from less complex cells**_
 - _Choosing appropriate thresholds for filtering, so as to **keep high quality cells without removing biologically relevant cell types**_

_**Recommendations:**_
 
 - _Have a good idea of your expectations for the **cell types to be present** prior to performing the QC. For instance, do you expect to have low complexity cells or cells with higher levels of mitochondrial expression in your sample? If so, then we need to account for this biology when filtering_
 - _If in doubt, be **less stringent with your filtering**. While having junk "cells" in the clustering can make it more difficult, if we see a cluster with high mitochondrial expression and no markers of interest, we can remove this cluster from the analysis downstream._

## Obtaining quality metrics for assessment

Oftentimes it is useful to extract the data and explore it ourselves without the limitations imposed by a specific tool/package. Therefore, we take the count matrix output and derive metrics for exploration. We have code [available](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_quality_control_analysis.html) for deriving these metrics if interested.

### Creating count data object

Generally, all single-cell RNA-seq datasets, regardless of technology or pipeline, will contain **three files**:

1. a file with the **gene IDs**, representing all genes quantified
2. a file with the **cell IDs**, representing all cells quantified
3. a **matrix of counts** per gene for every cell

You will need to use these files to properly give the matrix of counts the gene IDs as row names and cell IDs as column names.

Using this count matrix, we can create our metadata by calculating the different metrics to evaluate for quality control assessment. 

## QC metrics

We will be exploring multiple metrics used to determine which cells are of good quality. Below is an overview of the metrics they we often explore at the sample and cellular levels:

**Sample-level metrics include:**

- the distribution of number of reads per cell
- the number of cells per sample

These numbers should match our expectations based on what was asked for during sequencing. These metrics could help point to possible issues with specific samples.

For cellular quality control, we need to remove those "cells" that are really just free floating RNA that was encapsulated, as well as, dead or dying cells and low complexity cells like RBCs. To determine whether a cell is poor quality, there are several cellular-level metrics that we can explore.

**Cellular-level metrics include:**

- UMI counts per cell
- Genes detected per cell
- UMIs vs. genes detected
- Mitochondrial counts ratio
- Novelty

We will explore each of these metrics in more detail and discuss expectations and/or implications.

## Reads per cell

Generally, with this metric you hope to see all of the samples with peaks in the distributions to be in relatively the same location between 10,000 and 100,000 reads per cell. 

<p align="center">
<img src="../img/sc_qc_reads.png" width="800">
</p>
 
## Cell counts

The cell counts are determined by the number of unique cellular barcodes detected. **You expect the number of unique cellular barcodes to be around the number of sequenced cells or greater, sometimes much greater.** In single-cell protocols using hydrogels, like inDrops, some hydrogels may have more than one cellular barcode (see details in note below). After we remove the low quality cells by filtering, *we will expect the number of cells to be at or below the number of sequenced cells*. 

The capture efficiency of the method will determine how many of the sequenced cells are actually true high quality cells; inDrops tends to have a higher capture rate than 10X Genomics.

> **NOTE:** During the **inDrops** protocol, the cellular barcodes are present in the hydrogels, which are encapsulated in the droplets with a single cell and lysis/reaction mixture. Upon treatment of UV and cell lysis, all components mix together inside the droplet and reverse transcription proceeds, followed by droplet breakup and linear amplification for library preparation. While each hydrogel should have a single cellular barcode associated with it, occasionally a hydrogel can have more than one cellular barcode. We often see all possible combinations of cellular barcodes at a low level, leading to a higher number of cellular barcodes than cells.

<p align="center">
<img src="../img/sc_qc_cellcounts.png" width="500">
<\p>

## UMI counts (transcripts) per cell

The UMI counts per cell should generally be above 500, although usable, it's still low if between 500-1000 counts. If UMIs per cell is 500-1000 counts, then the cells probably should have been sequenced more deeply. 

<p align="center">
<img src="../img/sc_qc_umisPerCell.png" width="500">
<\p>
   
## Genes detected per cell

Seeing gene detection in the range of 500-5000 is normal for droplet-based analysis. Similar expectations for gene detection as for UMI detection, although may be a bit lower than UMIs.

<p align="center">
<img src="../img/sc_qc_genesDetected.png" width="500">
<\p>

## UMIs vs. genes detected

Poor quality cells are likely to have low genes and UMIs per cell. Therefore, a poor sample is likely to have cells in the lower left of the graph. Good cells should exhibit both higher number of genes per cell and higher numbers of UMIs. We also expect similar lines with similar slopes for all samples.

<p align="center">
<img src="../img/sc_qc_UMIsVsGenesDetected.png" width="500">
<\p>

## Mitochondrial counts ratio

This metric can identify whether there is a large amount of mitochondrial contamination from dead or dying cells. Poor quality samples for mitochondrial counts would have larger peaks above the 0.1-0.2 mitochondrial ratio mark, unless it is expected based on known cell types.

<p align="center">
<img src="../img/sc_qc_mitoRatio.png" width="500">
<\p>

## Novelty

We can see the samples where we sequenced each cell less have a higher overall novelty, that is because we have not started saturating the sequencing for any given gene for these samples. Outlier cells in these samples might be cells that have a less complex RNA species than other cells. Sometimes we can detect contamination with low complexity cell types like red blood cells via this metric. Generally, we expect the novelty score to be above 0.80.

<img src="../img/sc_qc_novelty.png" width="500">


## Filtering

Now that we have visualized the various metrics, we can decide on the thresholds to use to remove the low quality cells. Often the recommendations mentioned earlier are a rough guideline, but the specific experiment needs to inform the exact thresholds chosen. 

>**NOTE:** Non-droplet-based methods will have differing expectations for these metrics. For example, SMART-seq will return higher numbers of genes, but lower numbers of cells. SMART-seq methods will also not have UMIs present, so many of the metrics would not be used.

In addition to filtering the cells, we also filter to **remove lowly expressed genes**. Generally, we only keep genes that are expressed (counts > 0) in at least 10 cells.

After performing the filtering, it's recommended to look back over the metrics to make sure that your data matches your expectations and is good for downstream analysis. 

We will use this filtered count matrix for clustering our cells downstream, so it's important to remove any junk prior to this process.

[Click here for next lesson](https://github.com/marypiper/WIB_scRNA-seq/blob/master/lessons/SC_clustering_analysis.md)

---

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

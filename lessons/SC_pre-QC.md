---
title: "Single-cell RNA-seq: Generation of count matrix"
author: "Mary Piper, Lorena Pantano, Meeta Mistry, Radhika Khetani, Rory Kirchner"
date: Monday, April 22nd, 2019
---

Approximate time: 30 minutes

## Learning Objectives:

* Understand the general steps leading to generation of the count matrix

# Single-cell RNA-seq: raw sequencing data to counts

Single-cell RNA-seq (scRNA-seq) is an exciting and cutting-edge method for analyzing differences in cellular gene expression, particularly for tissue heterogeneity analyses, lineage tracing, and cell population dynamics. 

<img src="../img/sc_analyses.png" width="900">

The complexity of scRNA-seq data, which is generally characterized as:

- a **large volume of data:** expression data from thousands of cells
- a **low depth of sequencing per cell:** results in a large number of genes without any corresponding reads (zero inflation)

These characteristics make the **analysis of the data more involved** than bulk RNA-seq. In addition, the analyses can vary depending whether the goal is marker identification, lineage tracing, or some other custom analysis. Therefore, tools specific for scRNA-seq and the different methods of library preparation are needed. 

## Single-cell RNA-seq data

Depending on the library preparation method used, the RNA sequences, or reads/tags, will be derived either from the 3' ends (or 5' ends) of the transcripts (10X Genomics, CEL-seq2, Drop-seq, inDrop) or from full-length transcripts (Smart-seq). 

https://www.nature.com/articles/nri.2017.76

The choice of method involves the biological question of interest. The following advantages are listed below for the methods:

- **3' or 5'-end sequencing:** 
	- More accurate quantification through use of unique molecular identifiers distinguishing biological duplicates from amplification (PCR) duplicates
	- Larger number of cells sequenced allows better identity of cell type populations

- **Full length sequencing:**
	- Detection of isoform-level differences in expression
	- Identification of allele-specific differences in expression
	- Deeper sequencing of a smaller number of cells  

> **NOTE:** For the 3' sequencing method, reads originating from different molecules of the same transcript would likely have the same sequence (3' end). However, the PCR step during library preparation could also generate read duplicates. To determine whether a read is a biological or technical duplicate, these methods use unique molecular identifiers, or UMIs. Reads with different UMIs mapping to the same transcript were derived from different molecules and are biological duplicates, while reads with the same UMI originated from the same molecule and are technical duplicates.

## Single-cell RNA-seq workflow

The analysis workflow for scRNA-seq is generally similar for the differing scRNA-seq methods, but some specifics, such as the parsing of the UMIs, cell IDs, and sample IDs, will differ between them. For example, below is a schematic of the inDrop sequence reads:

<p align="center">
<img src="../img/sc_seq_method.png" width="600">
</p>

*Image credit: [Sarah Boswell](https://scholar.harvard.edu/saboswell), Director of the Single Cell Sequencing Core at HMS*

While the 10X sequence reads have the UMI and barcodes placed differently:

<p align="center">
<img src="../img/10_seq_method.png" width="600">
</p>

*Image credit: [Sarah Boswell](https://scholar.harvard.edu/saboswell), Director of the Single Cell Sequencing Core at HMS*

The scRNA-seq method will determine the how to parse the barcodes and UMIs from the sequencing reads. So, although a few of the specific steps will slightly differ, the overall workflow will generally follow the same steps regardless of method. The general workflow is shown below:

<img src="../img/sc_workflow.png" width="800">

The steps of the workflow are:

- **Generation of the count matrix (method-specific steps):** formating reads, demultiplexing samples, mapping and quantification
- **Quality control of the raw counts:** filtering of poor quality cells 
- **Clustering of filtered counts:** clustering cells based on similarities in transcriptional activity (cell types = different clusters)
- **Marker identification:** identifying gene markers for each cluster
- **Optional downstream steps**

## Generation of count matrix

We are going to start by discussing the first part of this workflow, which is generating the count matrix from the raw sequencing data. We will focus on the 3' end sequencing used by droplet-based methods, such as inDrop, 10X Genomics, and Drop-seq.

<p align="center">
<img src="../img/sc_gen_matrix_workflow.png" width="300">
</p>

After sequencing, the sequencing facility will either output the raw sequencing data as **BCL or FASTQ format**. If the reads are in BCL format, then we will need to convert to FASTQ format. There is a useful command-line tool called `bcl2fastq` that can easily perform this conversion. 

> **NOTE:** We do not demultiplex at this step in the workflow. You may have sequenced 6 samples, but the reads for all samples may be present all in the same BCL or FASTQ file.

The generation of the count matrix from the raw sequencing data will go through similar steps for many of the scRNA-seq methods. 

<img src="../img/sc_pre-QC_workflow.png" width="800">

[**umis**](https://github.com/vals/umis) and [**zUMIs**](https://github.com/sdparekh/zUMIs) are command-line tools that estimate expression of scRNA-seq data for which the 3' ends of transcripts were sequenced. Both tools incorporate collapsing of unique molecular identifiers to
correct for amplification bias. The steps in this process include the following:

 1. Formatting reads and filtering noisy cellular barcodes
 2. Demultiplexing the samples
 3. Pseudo-mapping to cDNAs
 4. Collapsing unique molecular identifiers and quantification of reads

## 1. Formatting reads and filtering noisy cellular barcodes

The FASTQ files can then be used to parse out the cell barcodes, UMIs, and sample barcodes. For droplet-based methods, many of the cellular barcodes will match a low number of reads (< 1000 reads) due to:

- encapsulation of free floating RNA from dying cells
- uncomplex small cells (RBCs, etc.)
- cells that failed for some reason

These excess barcodes need to be filtered out of the sequence data prior to read alignment. To do this filtering the 'cellular barcode' and the 'molecular barcode' is extracted and saved for each cell. For example, if using 'umis' tools, the information is added to the header line for each read, with the following format:

    @HWI-ST808:130:H0B8YADXX:1:1101:2088:2222:CELL_GGTCCA:UMI_CCCT
    AGGAAGATGGAGGAGAGAAGGCGGTGAAAGAGACCTGTAAAAAGCCACCGN
    +
    @@@DDBD>=AFCF+<CAFHDECII:DGGGHGIGGIIIEHGIIIGIIDHII#

Known cellular barcodes used in the library preparation method should be known, and unknown
barcodes would be dropped, however, allowing for an acceptable number of mismatches to the known cellular barcodes.

## 2. Demultiplexing sample reads

The next step of the process is to demultiplex the samples, if sequencing more than a single sample. This is the one step of this process not handled by the 'umis' tools, but is accomplished by 'zUMIs'. We would need to parse the reads to determine the sample barcode associated with each cell.

## 3. Pseudo-mapping to cDNAs

'This is done by pseudo-aligners, either Kallisto or RapMap. The SAM (or BAM) file output
from these tools need to be saved.'

## 4. Counting unique molecular identifiers

'The final step is to infer which cDNA was the origin of the tag a UMI was
attached to. We use the pseudo-alignments to the cDNAs, and consider a tag
assigned to a cDNA as a partial _evidence_ for a (cDNA, UMI) pairing. For
actual counting, we **only count unique UMIs** for (gene, UMI) pairings with
sufficient evidence.'

At this point of the workflow, the duplicate UMIs will be collapsed for the counting of the identifiers.

<img src="../img/sc_collapsing_umis.png" width="400">

Now we have our count matrix containing the counts per gene for each cell, which we can use to explore our data for quality information.

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

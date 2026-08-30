# RNA-seq Analysis of Drosophila melanogaster

## Overview

This repository contains the analysis workflow for an RNA-seq experiment
comparing untreated and RNAi-treated samples of *Drosophila melanogaster*.

The analysis was performed from raw sequencing reads through quality control,
adapter trimming, genome alignment, BAM processing, read counting, and
differential expression analysis.

The main biological comparison is:

RNAi vs Untreated

---

## Experimental Design

Six RNA-seq samples were analyzed.

| Sample | Condition |
|--------|-----------|
| SRR031708 | Untreated |
| SRR031712 | Untreated |
| SRR031713 | Untreated |
| SRR031720 | RNAi |
| SRR031721 | RNAi |
| SRR031722 | RNAi |

There are three biological samples in each condition.

---

## Reference Genome

Reference annotation:

Drosophila melanogaster BDGP5.25.62

The corresponding HISAT2 index and GTF annotation were used for
alignment and read counting.

---

## Analysis Pipeline

The workflow is organized into numbered scripts:

1. `01_fastqc.sh`  
   Initial quality control of raw FASTQ files.

2. `02_trimming.sh`  
   Adapter trimming and minimum read-length filtering.

3. `03_fastqc_trimmed.sh`  
   Quality control after trimming.

4. `04_hisat2_alignment.sh`  
   Alignment of trimmed reads to the Drosophila reference genome using HISAT2.

5. `05_bam_processing.sh`  
   Sorting, indexing and basic validation of BAM files using samtools.

6. `06_featurecounts.sh`  
   Assignment of aligned reads to genes using featureCounts.

7. `07_deseq2.R`  
   Differential expression analysis using DESeq2.

---

## Differential Expression Analysis

DESeq2 was used to compare:

RNAi vs Untreated

Low-count filtering:

- Count >= 10
- Present in at least 3 samples

Genes were considered differentially expressed using:

- adjusted p-value < 0.1
- absolute log2 fold change > 1

The final analysis identified:

- 258 differentially expressed genes
- 155 upregulated genes
- 103 downregulated genes

---

## Quality Control and Exploratory Analysis

The analysis includes:

- FastQC before trimming
- FastQC after trimming
- HISAT2 alignment statistics
- samtools BAM validation
- featureCounts assignment statistics
- PCA
- MA plot
- Volcano plot
- DEG heatmap

The PCA showed separation between the untreated and RNAi groups.

---

## Main Results

The complete DESeq2 results are available in:

`results/DEG/DESeq2_all_results.csv`

Annotated results:

`results/DEG/DESeq2_all_results_annotated.csv`

Final DEG table:

`results/DEG/DEG_final.csv`

Plots are available in:

`results/plots/`

---

## Software

Main software used in the analysis includes:

- FastQC
- cutadapt
- HISAT2 2.2.1
- samtools
- featureCounts 2.0.8
- R 4.3.3
- DESeq2 1.42.1
- ggplot2
- pheatmap

---

## Reproducibility

The numbered scripts in `scripts/` document the main computational steps
used for this analysis.

Large sequencing and alignment files are intentionally excluded from the
Git repository using `.gitignore`.

The repository therefore contains the analysis scripts, metadata,
counting results, differential-expression results and figures rather than
the original FASTQ and BAM files.

---

## Project Structure

```text
RNAseq/
├── scripts/
│   ├── 01_fastqc.sh
│   ├── 02_trimming.sh
│   ├── 03_fastqc_trimmed.sh
│   ├── 04_hisat2_alignment.sh
│   ├── 05_bam_processing.sh
│   ├── 06_featurecounts.sh
│   └── 07_deseq2.R
│
├── results/
│   ├── counts/
│   ├── DEG/
│   ├── plots/
│   └── metadata.csv
│
├── qc/
├── qc_trimmed/
├── reference/
├── .gitignore
└── README.md
```
## Notes

This repository documents the analysis workflow and the resulting files.

Raw sequencing data, trimmed FASTQ files, BAM files and local software
environments are not included in the Git repository because of their size.

The scripts are provided to make the analysis steps transparent and
reproducible.


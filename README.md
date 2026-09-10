# Bulkrnaseq-Analysis-Ulcerative-Colitis-vs-Healthy-Ileum-
Bulk RNA-seq differential expression analysis of Ulcerative Colitis vs Healthy ileum samples using limma-voom, with batch correction, DEG classification, and Reactome pathway enrichment (GSEA).

# Overview
This repository contains a complete Bulk RNA-seq analysis workflow for comparing gene expression profiles between Ulcerative Colitis (UC) and Healthy ileum samples.
The workflow starts with raw sequencing data and covers:
•	SRA data download
•	FASTQ conversion
•	FASTQ compression
•	Sequence quality control using FastQC
•	Multi-sample QC using MultiQC
•	Transcript quantification using Salmon
•	Generation of a raw count matrix
•	Metadata preparation
•	Protein-coding gene filtering
•	Batch-effect correction using limma
•	Principal Component Analysis (PCA)
•	Differential expression analysis using limma
•	DEG classification
•	Volcano plot visualization
•	Reactome pathway enrichment using GSEA

# Workflow
SRA accession IDs
        ↓
SRA Download
        ↓
FASTQ Conversion
        ↓
FASTQ Compression
        ↓
FastQC
        ↓
MultiQC
        ↓
Salmon Index
        ↓
Salmon Quantification
        ↓
quant.sf files
        ↓
Raw Count Matrix
        +
Metadata
        ↓
Protein-coding Gene Filtering
        ↓
Batch Correction
        ↓
PCA
        ↓
Differential Expression Analysis
        ↓
DEG Classification
        ↓
Volcano Plot
        ↓
Reactome GSEA
        ↓
Pathway Visualization

# 1. System and Software Requirements
The Linux portion of the workflow requires:
•	Ubuntu/Linux
•	Python 3
•	SRA Toolkit
•	pigz
•	FastQC
•	MultiQC
•	Salmon
The R analysis requires packages including:
•	tidyverse
•	limma
•	edgeR
•	clusterProfiler
•	org.Hs.eg.db
•	enrichplot
•	ggplot2
•	tibble
•	readr
•	ggrepel
•	factoextra

# 2. SRA Data Download
The sequencing datasets were obtained using SRA accession IDs.
Example accession IDs:
SRR23955797
SRR23955798
SRR23955799
Download an SRA dataset using:
prefetch SRR23955797
was performed for all samples.

# 3. FASTQ Conversion
The downloaded SRA files were converted into paired-end FASTQ files using fasterq-dump.
fasterq-dump SRR23955797.sra \
  --split-files \
  --threads 2 \
  --temp tmp \
  -O fastq
This generates paired FASTQ files such as:
SRR23955797_1.fastq
SRR23955797_2.fastq
The same procedure was repeated for each sample.

# 4. FASTQ Compression
FASTQ files were compressed using pigz to reduce storage requirements.
pigz -p 2 fastq/SRR23955797*.fastq
This produces:
SRR23955797_1.fastq.gz
SRR23955797_2.fastq.gz
Compression was performed for all samples.

# 5. Quality Control with FastQC
FastQC was used to assess the quality of the sequencing reads.
Example:
fastqc SRR23955797_1.fastq.gz \
  -o fastqc_results/
FastQC reports were generated for the sequencing samples.
Important QC metrics include:
•	Per-base sequence quality
•	Sequence quality scores
•	GC content
•	Sequence duplication
•	Adapter content
•	Overrepresented sequences

# 6. MultiQC
MultiQC was used to combine individual FastQC reports into a single summary report.
multiqc fastqc_results/ \
  -o multiqc_results/
The MultiQC report provides an overall view of sequencing quality across samples.

# 7. Salmon Transcript Quantification
Salmon was used for transcript-level quantification.
Build Salmon Index
A transcriptome FASTA file was used to generate the Salmon index:
salmon index \
  -t transcriptome.fa \
  -i salmon_index \
  -k 31
Quantify Samples
Each paired-end sample was quantified individually.
Example:
salmon quant \
  -i salmon_index \
  -l A \
  -1 SRR23955797_1.fastq.gz \
  -2 SRR23955797_2.fastq.gz \
  -p 12 \
  --validateMappings \
  -o SRR23955797_salmon_out
The same command was repeated for all samples.
Each Salmon output directory contains a quant.sf file containing transcript-level abundance estimates.
Example:
SRR23955797_salmon_out/
└── quant.sf

# 8. Generation of Raw Count Matrix
The individual Salmon quantification results were combined to generate a raw count matrix.
The resulting matrix contains:
•	Genes/transcripts as rows
•	Samples as columns
•	Expression/count values as entries
The main count matrix used for downstream analysis was:
ileum_raw_counts.csv

# 9. Metadata
A metadata file was prepared to describe the biological condition and batch information for each sample.
metadata.csv
The metadata contains information required for the downstream statistical analysis, including:
•	Sample ID
•	Condition
•	Batch
The sample names in the metadata were matched with the columns of the count matrix before analysis.

# 10. Protein-Coding Gene Filtering
Only protein-coding genes were retained for downstream analysis.
The list of protein-coding genes was stored in:
pc_gene.txt
The count matrix was filtered using this gene list.
protein_coding_genes <- read_lines("bulk/pc_gene.txt")

merged <- merged[
  rownames(merged) %in% protein_coding_genes,
]
Sample names were then checked to ensure that the count matrix and metadata were correctly matched.

# 11. Missing Gene / Sample QC
Missing or non-finite values were handled before downstream analysis.
merged[is.na(merged)] <- 0
merged[!is.finite(as.matrix(merged))] <- 0
The proportion of zero-count genes was calculated for each sample.
sample_missing_pct <-
  colSums(merged == 0) / nrow(merged)
The results were saved as:
Percentage_of_missing_genes_per_samples.csv
A histogram was also generated to visualize the distribution of missing genes across samples.

# 12. Batch Correction
Batch effects were assessed and corrected using the limma package.
The design matrix was constructed using biological condition:
mod <- model.matrix(~ Condition, data = metadata)
Batch correction was performed using:
corrected_data_limma_bc <-
  removeBatchEffect(
    as.matrix(merged),
    batch = metadata$Batch,
    design = mod
  )
Two datasets were therefore available for comparison:
•	Before batch correction
•	After batch correction
The corrected expression matrix was saved as:
limma_merged.csv

# 13. PCA Before and After Batch Correction
Principal Component Analysis (PCA) was performed to evaluate sample-level structure and the effect of batch correction.
PCA was visualized according to:
•	Biological condition
•	Batch
The following plots were generated:
PCA_before_Condition.png
PCA_after_Condition.png

PCA_before_Batch.png
PCA_after_Batch.png
These plots were used to evaluate whether batch-related separation was reduced after correction while retaining biological differences.

# 14. Differential Expression Analysis
Differential expression analysis was performed using limma.
The corrected expression matrix was loaded:
ileum_data <-
  read.csv(
    "bulk/output/limma_merged.csv",
    row.names = 1
  )
The experimental groups were defined as:
group <- factor(metadata$Condition)

design <- model.matrix(~0 + group)

colnames(design) <- levels(group)
The linear model was fitted:
fit <- lmFit(ileum_data, design)
The comparison was defined as:
UC_No vs Healthy
using:
contrast.matrix <-
  makeContrasts(
    DiseaseUntreated_vs_Healthy =
      UC_No - Healthy,
    levels = design
  )
The model was then fitted and moderated using empirical Bayes statistics:
fit2 <- contrasts.fit(
  fit,
  contrast.matrix
)

fit2 <- eBayes(fit2)
The complete differential expression results were obtained using:
top_genes <-
  topTable(
    fit2,
    adjust = "fdr",
    number = Inf
  )
The results were saved as:
DEGsList.csv

# 15. DEG Classification
Genes with:
Adjusted P value < 0.05
were considered significant DEGs.
The significant genes were classified into:
•	Upregulated genes
•	Downregulated genes
The total number of DEGs and the numbers of upregulated and downregulated genes were calculated.

# 16. PCA of UC vs Healthy Samples
A PCA plot was generated to visualize separation between the two biological conditions.
The resulting figure was saved as:
contrastPCA.png

# 17. Volcano Plot
A volcano plot was generated to visualize differential expression.
The plot displays:
•	Log2 fold change on the X-axis
•	−log10 adjusted P value on the Y-axis
Genes were categorized as:
•	Not significant
•	Significant
•	Reported in Papers
Selected genes of interest were highlighted and labelled.
The highlighted genes included:
REG1A
FDFT1
SIK1
DMBT1
GGACT
DUSP21
The final figure was saved as:
VolcanoPlot.png

# 18. Reactome GSEA
Reactome pathway enrichment was performed using a ranked gene list.
Gene symbols were converted to Entrez IDs using:
bitr(
  rownames(res),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
The genes were ranked according to their log fold change:
ngenes <- res_mapped$logFC

names(ngenes) <- res_mapped$ENTREZID

ngenes <- sort(
  ngenes,
  decreasing = TRUE
)
Reactome Gene Set Enrichment Analysis was then performed:
enp_gsea <- gsePathway(
  geneList = ngenes,
  organism = "human",
  verbose = FALSE
)
The enrichment results were converted to readable gene identifiers and ranked according to adjusted P value and normalized enrichment score.

# 19. Reactome Pathway Visualization
The top 20 Reactome pathways were visualized using a dot plot.
The plot represents:
•	Normalized Enrichment Score (NES)
•	Adjusted P value
•	Gene set size
The results were saved as:
pathways.csv
and the visualization as:
pathwaysplot.png

# 20. Repository Structure
Bulk-RNAseq-UC-Ileum/
│
├── README.md
│
├── data/
│   ├── metadata.csv
│   ├── ileum_raw_counts.csv
│   └── pc_gene.txt
│
├── scripts/
│   └── BulkRNAseq_analysis.R
│
├── results/
│   ├── QC/
│   │   ├── Percentage_of_missing_genes_per_samples.csv
│   │   └── Histogram_percentage_of_missing_values_per_samples.png
│   │
│   ├── Batch_Correction/
│   │   ├── limma_merged.csv
│   │   ├── PCA_before_Condition.png
│   │   ├── PCA_after_Condition.png
│   │   ├── PCA_before_Batch.png
│   │   └── PCA_after_Batch.png
│   │
│   ├── Differential_Expression/
│   │   ├── DEGsList.csv
│   │   ├── contrastPCA.png
│   │   └── VolcanoPlot.png
│   │
│   └── Reactome_GSEA/
│       ├── pathways.csv
│       └── pathwaysplot.png
│
└── .gitignore

# 21. Key Tools
Tool	Purpose
SRA Toolkit	SRA download and FASTQ conversion
pigz	FASTQ compression
FastQC	Sequencing quality control
MultiQC	Aggregated QC reporting
Salmon	Transcript quantification
R	Statistical analysis
limma	Batch correction and differential expression
edgeR	RNA-seq analysis support
clusterProfiler	Gene ID conversion and enrichment analysis
ReactomePA	Reactome pathway GSEA
ggplot2	Data visualization

# 22. Key Outputs
The major outputs of this project include:
Raw count matrix
        ↓
Batch-corrected expression matrix
        ↓
PCA plots
        ↓
Differential expression results
        ↓
DEG classification
        ↓
Volcano plot
        ↓
Reactome GSEA results
        ↓
Top Reactome pathway visualization

# Project Summary
This project demonstrates a complete Bulk RNA-seq workflow for investigating transcriptional differences between Ulcerative Colitis and Healthy ileum samples.
Starting from sequencing data, the workflow performs quality control, transcript quantification, count matrix preparation, metadata integration, batch-effect correction, PCA, differential expression analysis, DEG classification, and Reactome pathway enrichment.
The analysis provides both gene-level differential expression results and pathway-level biological interpretation of transcriptional changes associated with Ulcerative Colitis.


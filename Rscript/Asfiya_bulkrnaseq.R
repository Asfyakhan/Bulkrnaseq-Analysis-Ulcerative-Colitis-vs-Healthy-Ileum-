# Bulk RNA-seq Analysis: Ulcerative Colitis vs Healthy (Ileum)

###Chunk 1: SETUP: Install required packages###

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages(c("tidyverse", "ggrepel", "factoextra", "readxl"))
BiocManager::install(c(
  "limma", "edgeR", "DESeq2", "clusterProfiler", "org.Hs.eg.db",
  "enrichplot", "ReactomePA", "DOSE", "rtracklayer"
))

setwd("D:/BulkRNAseq/Deseq")
dir.create("output", showWarnings = FALSE)

###Chunk 2: Load metadata and raw counts###

detect_and_read <- function(pattern) {
  f <- list.files(pattern = pattern, full.names = TRUE, ignore.case = TRUE)[1]
  ext <- tolower(tools::file_ext(f))
  if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(f)
  } else {
    read.csv(f, sep = "\t", check.names = FALSE)
  }
}

getwd()
list.files()

detect_and_read <- function(pattern) {
  matches <- list.files(pattern = pattern, ignore.case = TRUE, full.names = TRUE)
  
  if (length(matches) == 0) {
    stop("No file found matching pattern '", pattern, "' in: ", getwd(),
         "\nFiles actually here: ", paste(list.files(), collapse = ", "))
  }
  
  f <- matches[1]
  ext <- tolower(tools::file_ext(f))
  
  if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(f)
  } else {
    read.csv(f, sep = "\t", check.names = FALSE)
  }
}
metadata <- detect_and_read("^metadata")
merged   <- detect_and_read("^ileum_raw_counts")

dim(metadata)
dim(merged)
head(metadata)

#Align metadata column

library(tidyverse)

metadata <- metadata %>% rename(Batch = Dataset)
merged   <- merged[, colnames(merged) != ""]
rownames(merged) <- make.unique(rownames(merged))

dim(metadata)
dim(merged)
colnames(metadata)

###Chunk 3: Keep only protein-coding genes###

#Build protein-coding gene list from GTF

protein_coding_genes <- read_lines("pc_gene.txt")
length(protein_coding_genes)

list.files()
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("rtracklayer")
library(rtracklayer)
install.packages("RCurl")
library(rtracklayer)

gtf <- rtracklayer::import("pc_gene.txt.gtf")
gtf_df <- as.data.frame(gtf)

biotype_col <- intersect(
  c("gene_type", "gene_biotype"),
  colnames(gtf_df)
)[1]

protein_coding_genes <- gtf_df %>%
  filter(
    type == "gene",
    .data[[biotype_col]] == "protein_coding"
  ) %>%
  pull(gene_name) %>%
  unique()

length(protein_coding_genes)

#Filter Protein-Coding Genes
write_lines(protein_coding_genes, "pc_gene.txt")
merged <- merged[rownames(merged) %in% protein_coding_genes, ]
dim(merged)

# Check gene IDs in count matrix
head(rownames(merged), 10)

# Check protein-coding gene IDs
head(protein_coding_genes, 10)

# Check how many IDs match
sum(rownames(merged) %in% protein_coding_genes)

# Reload the original count matrix
merged <- detect_and_read("^ileum_raw_counts")

# Set first column as gene names
rownames(merged) <- merged[[1]]
merged[[1]] <- NULL

# Remove empty column names if present
merged <- merged[, colnames(merged) != ""]

# Make gene names unique
rownames(merged) <- make.unique(rownames(merged))

# Check dimensions
dim(merged)

# Check first 10 gene IDs
head(rownames(merged), 10)

# Check matching genes
sum(rownames(merged) %in% protein_coding_genes)

merged <- merged[rownames(merged) %in% protein_coding_genes, ]

dim(merged)

###Chunk 4: Match metadata to counts & clean groups###

#Match Metadata and Sample Order

intersect(colnames(merged), metadata$Sample)

metadata <- metadata[metadata$Sample %in% colnames(merged), ]

metadata <- metadata[match(colnames(merged), metadata$Sample), ]

all(metadata$Sample == colnames(merged))

#Check common sample names
intersect(colnames(merged), metadata$Sample)

#Keep only metadata samples present in count matrix
metadata <- metadata[metadata$Sample %in% colnames(merged), ]

#Reorder metadata to exactly match count matrix columns
metadata <- metadata[match(colnames(merged), metadata$Sample), ]

#Final check - must return TRUE
all(metadata$Sample == colnames(merged))

###Chunk 5: QC missing-gene histogram & save cleaned raw counts###

merged[is.na(merged)] <- 0
merged[!is.finite(as.matrix(merged))] <- 0

sample_missing_pct <- colSums(merged == 0) / nrow(merged)
missing_df <- data.frame(Sample = colnames(merged),
                         Missing_Percentage = sample_missing_pct[colnames(merged)])
write.csv(missing_df, "Percentage_of_missing_genes_per_samples.csv", row.names = FALSE)
write.csv(missing_df, "Percentage_of_missing_genes_per_samples.csv", row.names = FALSE)

list.files()
head(missing_df)
summary(missing_df$Missing_Percentage)

sample_missing_pct <- colSums(merged == 0) / nrow(merged)

missing_df <- data.frame(
  Sample = colnames(merged),
  Missing_Percentage = sample_missing_pct * 100
)

print(missing_df)

write.csv(
  missing_df,
  "Percentage_of_missing_genes_per_samples.csv",
  row.names = FALSE
)

#Save histogram

if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}

#QC HISTOGRAM

png(
  "output/Histogram_percentage_of_missing_values_per_samples.png",
  width = 10,
  height = 8,
  units = "in",
  res = 300
)

hist(
  missing_df$Missing_Percentage,
  main = "Missing Genes per Sample",
  xlab = "Percentage of Genes with Zero Counts",
  ylab = "Number of Samples"
)

dev.off()

list.files("output")


#SAVE PROCESSED COUNTS


write.csv(
  merged,
  "output/raw_merged.csv",
  row.names = TRUE
)

dim(merged)

###Chunk 6: log-CPM & PCA before batch correction###

#Prepare metadata and check batch structure

#PREPARE METADATA

setwd("D:/BulkRNAseq/Deseq")
list.files()

getwd()
install.packages("readxl")
library(readxl)

metadata <- read_excel("D:/BulkRNAseq/metadata_csv.xlsx")


metadata$Batch <- as.factor(metadata$Dataset)

metadata$Condition <- as.factor(metadata$Condition)

metadata <- metadata[
  match(colnames(merged), metadata$Sample),
]

stopifnot(all(metadata$Sample == colnames(merged)))

head(metadata)

ls()

exists("merged")
exists("metadata")
exists("missing_df")
exists("protein_coding_genes")

setwd("D:/BulkRNASeq/Deseq")
getwd()
library(readxl)
metadata <- read_excel(
  "D:/BulkRNASeq/Deseq/metadata.csv.xlsx"
)

head(metadata)
colnames(metadata)
dim(metadata)
merged <- read.csv(
  "D:/BulkRNASeq/Deseq/ileum_raw_counts.csv",
  check.names = FALSE,
  row.names = 1
)
dim(merged)
head(merged[, 1:min(5, ncol(merged))])
rm(merged)

count_test <- read.csv(
  "D:/BulkRNASeq/Deseq/ileum_raw_counts.csv",
  check.names = FALSE
)

dim(count_test)
colnames(count_test)
readLines(
  "D:/BulkRNASeq/Deseq/ileum_raw_counts.csv",
  n = 3
)
readLines(
  "D:/BulkRNASeq/Deseq/ileum_raw_counts.csv",
  n = 3
)
count_test2 <- read.csv(
  "D:/BulkRNASeq/Deseq/ileum_raw_counts.csv",
  sep = ";",
  check.names = FALSE
)

dim(count_test2)
count_test3 <- read.delim(
  "D:/BulkRNASeq/Deseq/ileum_raw_counts.csv",
  check.names = FALSE
)

dim(count_test3)

colnames(count_test3)[1:10]

head(count_test3[, 1:6])

str(count_test3[, 1:6])

head(count_test3[[1]])

#Save filtered raw count matrix

write.csv(
  merged,
  "output/raw_merged.csv"
)

count_test3 <- read.delim(
  "D:/BulkRNASeq/Deseq/ileum_raw_counts.csv",
  check.names = FALSE
)

merged <- count_test3

rownames(merged) <- make.unique(
  as.character(merged[[1]])
)

merged <- merged[, -1, drop = FALSE]
dim(merged)
length(protein_coding_genes)
protein_coding_genes <- readLines(
  "D:/BulkRNASeq/Deseq/pc_gene.txt"
)
length(protein_coding_genes)
exists("merged")
dim(merged)

#Save filtered raw count matrix

write.csv(
  merged,
  "output/raw_merged.csv"
)
file.exists("output/raw_merged.csv")

# Make sure metadata order matches merged columns
metadata <- metadata[
  match(colnames(merged), metadata$Sample),
  ,
  drop = FALSE
]

# Stop if sample order is incorrect
stopifnot(
  all(metadata$Sample == colnames(merged))
)

# Convert Batch to factor
metadata$Batch <- as.factor(metadata$Batch)

# Convert Condition to factor
metadata$Condition <- as.factor(metadata$Condition)

# Create design matrix using Condition
mod <- model.matrix(
  ~ Condition,
  data = metadata
)

colnames(metadata)
names(metadata)

"Dataset" %in% colnames(metadata)

"Condition" %in% colnames(metadata)

"Sample" %in% colnames(metadata)

metadata$Batch <- as.factor(metadata$Dataset)

metadata$Condition <- as.factor(metadata$Condition)

table(metadata$Batch)

table(metadata$Condition)

table(metadata$Batch, metadata$Condition)


exclude <- c(
  "GSE183620",
  "GSE193677"
)

# Drop unused batch levels
metadata$Batch <- droplevels(metadata$Batch)

# Remove excluded batches
metadata <- metadata[
  !metadata$Batch %in% exclude,
  ,
  drop = FALSE
]

# Keep only corresponding samples
merged <- merged[
  ,
  colnames(merged) %in% metadata$Sample,
  drop = FALSE
]

# Reorder metadata to match count matrix
metadata <- metadata[
  match(colnames(merged), metadata$Sample),
  ,
  drop = FALSE
]

# Convert Batch to factor again
metadata$Batch <- as.factor(metadata$Batch)

# Convert Condition to factor
metadata$Condition <- as.factor(metadata$Condition)

# Check sample order
all(metadata$Sample == colnames(merged))

table(metadata$Batch)

table(metadata$Condition)

table(metadata$Batch, metadata$Condition)

#Final batch correction

library(limma)

# Create new design matrix
mod <- model.matrix(
  ~ Condition,
  data = metadata
)

# Perform final batch correction
corrected_data_limma_bc <- removeBatchEffect(
  as.matrix(merged),
  batch = metadata$Batch,
  design = mod
)

# Check dimensions
dim(corrected_data_limma_bc)

nrow(metadata)
ncol(merged)

unique(metadata$Dataset)
sum(metadata$Dataset %in% c("GSE183620", "GSE193677"))

keep <- metadata$Condition != "UC_Uninflamed"

metadata <- metadata[keep, , drop = FALSE]

merged <- merged[
  ,
  metadata$Sample,
  drop = FALSE
]
metadata <- metadata[
  match(colnames(merged), metadata$Sample),
  ,
  drop = FALSE
]
table(metadata$Condition)
metadata$Condition <- droplevels(
  as.factor(metadata$Condition)
)
table(
  metadata$Batch,
  metadata$Condition
)

mod <- model.matrix(
  ~ Condition,
  data = metadata
)

mod

metadata$Batch <- droplevels(
  as.factor(metadata$Batch)
)

library(limma)

corrected_data_limma_bc <- removeBatchEffect(
  as.matrix(merged),
  batch = metadata$Batch,
  design = mod
)

setwd("D:/BulkRNAseq/Deseq")
getwd()
list.files()
library(tidyverse)
metadata <- read.csv("metadata.csv.xlsx")
merged <- read.csv("ileum_raw_counts.csv", check.names = FALSE, row.names = 1)

metadata
merged
dim(metadata)
dim(merged)
library(readxl)

metadata <- read_excel("metadata.csv.xlsx")
metadata <- as.data.frame(metadata)

merged <- read.delim("ileum_raw_counts.csv", check.names = FALSE, row.names = 1)

metadata
merged
dim(metadata)
dim(merged)

metadata$Batch <- as.factor(metadata$Dataset)

metadata <- metadata[match(colnames(merged), metadata$Sample), ]
stopifnot(all(metadata$Sample == colnames(merged)))

metadata$Condition <- as.factor(metadata$Condition)
mod <- model.matrix(~ Condition, data = metadata)

library(limma)
corrected_data_limma_bc <- removeBatchEffect(as.matrix(merged),
                                             batch = metadata$Batch,
                                             design = mod)
levels(metadata$Batch)
table(metadata$Batch, metadata$Condition)

exclude <- c("GSE183620", "GSE193677")

metadata <- metadata[metadata$Condition %in% c("Healthy", "UC"), ]
metadata$Condition <- droplevels(as.factor(metadata$Condition))
metadata$Batch <- droplevels(as.factor(metadata$Batch))

merged <- merged[, colnames(merged) %in% metadata$Sample]
metadata <- metadata[match(colnames(merged), metadata$Sample), ]

stopifnot(all(metadata$Sample == colnames(merged)))

table(metadata$Batch, metadata$Condition)

# 3. Rebuild design matrix and re-run batch correction
mod <- model.matrix(~ Condition, data = metadata)

library(limma)
corrected_data_limma_bc <- removeBatchEffect(as.matrix(merged),
                                             batch = metadata$Batch,
                                             design = mod)

sum(is.na(corrected_data_limma_bc))

write.csv(corrected_data_limma_bc, "output/limma_merged.csv")

###Chunk 6 — PCA before vs after batch correction###

library(ggplot2)

merged <- merged[apply(merged, 1, var) > 0, ]
pca_before <- prcomp(t(merged), scale. = TRUE)
pca_before_df <- data.frame(pca_before$x, Condition = metadata$Condition, Batch = metadata$Batch)
var_explained_before <- round(100 * (pca_before$sdev^2 / sum(pca_before$sdev^2)), 2)


###Chunk 7: Batch correction & PCA after###

corrected_data_limma_bc <- corrected_data_limma_bc[apply(corrected_data_limma_bc, 1, var) > 0, ]
pca_after <- prcomp(t(corrected_data_limma_bc), scale. = TRUE)
pca_after_df <- data.frame(pca_after$x, Condition = metadata$Condition, Batch = metadata$Batch)
var_explained_after <- round(100 * (pca_after$sdev^2 / sum(pca_after$sdev^2)), 2)

plot_pca <- function(pca_df, var_exp, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"), color = "Condition") +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
}

plot_pcab <- function(pca_df, var_exp, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"), color = "Batch") +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
}

p3 <- plot_pca(pca_before_df, var_explained_before, "PCA Global Before - by Condition")
p4 <- plot_pca(pca_after_df, var_explained_after, "PCA Global Final - by Condition")
p2 <- plot_pcab(pca_before_df, var_explained_before, "PCA Global Before - by Batch")
p1 <- plot_pcab(pca_after_df, var_explained_after, "PCA Global Final - by Batch")

ggsave("output/PCA_before_Condition.png", p3, width = 8, height = 6, dpi = 300)
ggsave("output/PCA_after_Condition.png", p4, width = 8, height = 6, dpi = 300)
ggsave("output/PCA_before_Batch.png", p2, width = 8, height = 6, dpi = 300)
ggsave("output/PCA_after_Batch.png", p1, width = 8, height = 6, dpi = 300)

# Flag PCA outliers using PC1/PC2 from the corrected data (pca_after)

outlier_scores <- data.frame(
  Sample = rownames(pca_after$x),
  PC1 = pca_after$x[, "PC1"],
  PC2 = pca_after$x[, "PC2"]
)
outlier_scores$PC1_z <- scale(outlier_scores$PC1)
outlier_scores$PC2_z <- scale(outlier_scores$PC2)

# Samples more than 3 SD away on either axis

pca_outliers <- outlier_scores[abs(outlier_scores$PC1_z) > 3 | abs(outlier_scores$PC2_z) > 3, ]
pca_outliers

# Cross-check against your missing-genes QC metric from earlier

missing_df <- read.csv("output/Percentage_of_missing_genes_per_samples.csv", row.names = 1)
head(missing_df)
merge(pca_outliers, missing_df, by = "Sample")

write.csv(missing_df, "Percentage_of_missing_genes_per_samples.csv")

missing_df <- read.csv("Percentage_of_missing_genes_per_samples.csv", row.names = 1)
head(missing_df)

merge(pca_outliers, missing_df, by = "Sample")

colnames(pca_outliers)
colnames(missing_df)

raw_check <- read.csv("Percentage_of_missing_genes_per_samples.csv")
head(raw_check)
colnames(raw_check)

missing_df <- read.csv("Percentage_of_missing_genes_per_samples.csv")
colnames(missing_df)   
merge(pca_outliers, missing_df, by = "Sample")

outlier_ids <- c("GSM5035183", "GSM7114724", "GSM7114726", "GSM7114741")
metadata[metadata$Sample %in% outlier_ids, c("Sample", "Dataset", "Condition")]

library_sizes <- colSums(merged[, colnames(merged) %in% outlier_ids])
library_sizes
summary(colSums(merged))

library(edgeR)

# Convert raw counts to log2 CPM (normalizes for library size)

logcpm <- cpm(as.matrix(merged), log = TRUE, prior.count = 1)

pca_logcpm <- prcomp(t(logcpm), scale. = TRUE)
var_explained_logcpm <- round(100 * (pca_logcpm$sdev^2 / sum(pca_logcpm$sdev^2)), 2)

pca_logcpm_df <- data.frame(pca_logcpm$x,
                            Sample = colnames(logcpm),
                            Condition = metadata$Condition,
                            Batch = metadata$Batch)
library(ggplot2)
ggplot(pca_logcpm_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = paste0("PCA on log2-CPM (PC1: ", var_explained_logcpm[1],
                      "%, PC2: ", var_explained_logcpm[2], "%)"))

pca_logcpm_df[pca_logcpm_df$Sample %in% outlier_ids, c("Sample", "Batch", "PC1", "PC2")]

library(edgeR)
logcpm <- cpm(as.matrix(merged), log = TRUE, prior.count = 1)

metadata$Batch <- as.factor(metadata$Batch)
metadata <- metadata[match(colnames(logcpm), metadata$Sample), ]
stopifnot(all(metadata$Sample == colnames(logcpm)))

mod <- model.matrix(~ Condition, data = metadata)
library(limma)
corrected_data_limma_bc <- removeBatchEffect(logcpm,
                                             batch = metadata$Batch,
                                             design = mod)

sum(is.na(corrected_data_limma_bc))   
write.csv(corrected_data_limma_bc, "output/limma_merged.csv")

pca_before <- prcomp(t(logcpm), scale. = TRUE)
pca_before_df <- data.frame(pca_before$x, Condition = metadata$Condition, Batch = metadata$Batch)
var_explained_before <- round(100 * (pca_before$sdev^2 / sum(pca_before$sdev^2)), 2)

pca_after <- prcomp(t(corrected_data_limma_bc), scale. = TRUE)
pca_after_df <- data.frame(pca_after$x, Condition = metadata$Condition, Batch = metadata$Batch)
var_explained_after <- round(100 * (pca_after$sdev^2 / sum(pca_after$sdev^2)), 2)

library(ggplot2)

plot_pca <- function(pca_df, var_exp, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"), color = "Condition") +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
}

plot_pcab <- function(pca_df, var_exp, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"), color = "Batch") +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
}

p3 <- plot_pca(pca_before_df, var_explained_before, "PCA (log-CPM) Before - by Condition")
p4 <- plot_pca(pca_after_df, var_explained_after, "PCA (log-CPM) After - by Condition")
p2 <- plot_pcab(pca_before_df, var_explained_before, "PCA (log-CPM) Before - by Batch")
p1 <- plot_pcab(pca_after_df, var_explained_after, "PCA (log-CPM) After - by Batch")

ggsave("output/PCA_before_Condition.png", p3, width = 8, height = 6, dpi = 300)
ggsave("output/PCA_after_Condition.png", p4, width = 8, height = 6, dpi = 300)
ggsave("output/PCA_before_Batch.png", p2, width = 8, height = 6, dpi = 300)
ggsave("output/PCA_after_Batch.png", p1, width = 8, height = 6, dpi = 300)

#24/08/2026

setwd("D:/BulkRNAseq/Deseq")
getwd()
list.files()

if (!dir.exists("output")) dir.create("output")
library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(DESeq2)

library(readxl)

metadata <- read_excel("metadata.csv.xlsx", sheet = "data")
metadata <- as.data.frame(metadata)

merged <- read.delim("ileum_raw_counts.csv", check.names = FALSE, row.names = 1)

dim(metadata)
dim(merged)

rownames(merged) <- make.unique(rownames(merged))

protein_coding_file <- "pc_gene.txt"
protein_coding_genes <- read_lines(protein_coding_file)
merged <- merged[rownames(merged) %in% protein_coding_genes, ]

library(ggplot2)

plot_pca <- function(pca_df, var_exp, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"), color = "Condition") +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
}

plot_pcab <- function(pca_df, var_exp, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"), color = "Batch") +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
}

exists("pca_after_df")
exists("var_explained_after")

p1 <- plot_pcab(pca_after_df, var_explained_after, "PCA (log-CPM) After - by Batch")
ggsave("output/PCA_after_Batch.png", p1, width = 8, height = 6, dpi = 300)


setwd("D:/BulkRNAseq/Deseq")
list.files()

library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(ggplot2)

metadata <- read_excel("metadata.csv.xlsx", sheet = "data")   
metadata <- as.data.frame(metadata)

merged <- read.delim("ileum_raw_counts.csv", check.names = FALSE, row.names = 1)

rownames(merged) <- make.unique(rownames(merged))
protein_coding_genes <- read_lines("pc_gene.txt")
merged <- merged[rownames(merged) %in% protein_coding_genes, ]

metadata <- metadata[metadata$Sample %in% colnames(merged), ]
metadata <- metadata[match(colnames(merged), metadata$Sample), ]
stopifnot(all(metadata$Sample == colnames(merged)))

merged[is.na(merged)] <- 0
merged[!is.finite(as.matrix(merged))] <- 0

exclude <- c("GSE183620", "GSE193677")

metadata <- metadata[metadata$Condition %in% c("Healthy", "UC"), ]
metadata$Condition <- droplevels(as.factor(metadata$Condition))
metadata$Batch <- droplevels(as.factor(metadata$Dataset))

merged <- merged[, colnames(merged) %in% metadata$Sample]
metadata <- metadata[match(colnames(merged), metadata$Sample), ]
stopifnot(all(metadata$Sample == colnames(merged)))

exclude <- c("GSE183620", "GSE193677")

metadata <- metadata[metadata$Condition %in% c("Healthy", "UC"), ]
metadata$Condition <- droplevels(as.factor(metadata$Condition))
metadata$Batch <- droplevels(as.factor(metadata$Dataset))

merged <- merged[, colnames(merged) %in% metadata$Sample]
metadata <- metadata[match(colnames(merged), metadata$Sample), ]
stopifnot(all(metadata$Sample == colnames(merged)))

logcpm <- cpm(as.matrix(merged), log = TRUE, prior.count = 1)

mod <- model.matrix(~ Condition, data = metadata)
corrected_data_limma_bc <- removeBatchEffect(logcpm,
                                             batch = metadata$Batch,
                                             design = mod)
sum(is.na(corrected_data_limma_bc))  

pca_before <- prcomp(t(logcpm), scale. = TRUE)
pca_before_df <- data.frame(pca_before$x, Condition = metadata$Condition, Batch = metadata$Batch)
var_explained_before <- round(100 * (pca_before$sdev^2 / sum(pca_before$sdev^2)), 2)

pca_after <- prcomp(t(corrected_data_limma_bc), scale. = TRUE)
pca_after_df <- data.frame(pca_after$x, Condition = metadata$Condition, Batch = metadata$Batch)
var_explained_after <- round(100 * (pca_after$sdev^2 / sum(pca_after$sdev^2)), 2)

plot_pca <- function(pca_df, var_exp, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"), color = "Condition") +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
}

plot_pcab <- function(pca_df, var_exp, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"), color = "Batch") +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
}

p1 <- plot_pcab(pca_after_df, var_explained_after, "PCA (log-CPM) After - by Batch")
ggsave("output/PCA_after_Batch.png", p1, width = 8, height = 6, dpi = 300)

var_explained_after

dim(merged)      # genes x samples going into logcpm
dim(metadata)   # should match ncol(merged)
nrow(logcpm)     # genes used for PCA
ncol(logcpm)     # samples used for PCA

length(protein_coding_genes)
nrow(merged)   # should be close to length(protein_coding_genes), not 46105 (the unfiltered count)

table(metadata$Batch, metadata$Condition)

sum(apply(logcpm, 1, var) == 0)
sum(apply(corrected_data_limma_bc, 1, var) == 0)

###Chunk 8: DGEList, filter low counts, TMM normalize###

#differential expression analysis

library(tidyverse)
library(readxl)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
library(ggrepel)

write.csv(
  corrected_data_limma_bc,
  "output/logCPM_batch_corrected_limma.csv"
)

write.csv(
  metadata,
  "output/metadata_final_137_samples.csv",
  row.names = FALSE
)

write.csv(
  logcpm,
  "output/logCPM_before_batch_correction.csv"
)

lib_sizes <- colSums(merged)

summary(lib_sizes)

libsize_df <- data.frame(
  Sample = colnames(merged),
  Library_Size = lib_sizes,
  Condition = metadata$Condition,
  Batch = metadata$Batch
)

p_libsize <- ggplot(
  libsize_df,
  aes(x = Batch, y = Library_Size, fill = Condition)
) +
  geom_boxplot(alpha = 0.8) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Library Size Distribution by Batch",
    x = "Batch",
    y = "Total Raw Counts"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_libsize)

ggsave(
  "output/Library_Size_by_Batch.png",
  p_libsize,
  width = 9,
  height = 6,
  dpi = 300
)

dge <- DGEList(
  counts = merged,
  group = metadata$Condition
)
dge

#Filter lowly expressed genes

design <- model.matrix(
  ~ Batch + Condition,
  data = metadata
)

head(design)
dim(design)

keep <- filterByExpr(
  dge,
  design = design
)

table(keep)

dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

dim(dge_filtered)

#TMM normalization

dge_filtered <- calcNormFactors(
  dge_filtered,
  method = "TMM"
)

dge_filtered$samples


normalization_factors <- data.frame(
  Sample = rownames(dge_filtered$samples),
  dge_filtered$samples
)

write.csv(
  normalization_factors,
  "output/TMM_Normalization_Factors.csv",
  row.names = FALSE
)

#MDS plot after TMM normalization

png(
  "output/MDS_by_Condition.png",
  width = 2400,
  height = 1800,
  res = 300
)

plotMDS(
  dge_filtered,
  labels = metadata$Condition,
  col = ifelse(
    metadata$Condition == "Healthy",
    "tomato",
    "turquoise4"
  )
)

legend(
  "topright",
  legend = levels(metadata$Condition),
  col = c("tomato", "turquoise4"),
  pch = 16
)

dev.off()

#Recreate the design matrix carefully

metadata$Condition <- relevel(
  factor(metadata$Condition),
  ref = "Healthy"
)

metadata$Batch <- factor(metadata$Batch)

design <- model.matrix(
  ~ Batch + Condition,
  data = metadata
)

colnames(design)

qr(design)$rank
ncol(design)

###Chunk 9: voom transformation###

v <- voom(
  dge_filtered,
  design,
  plot = TRUE
)

png(
  "output/Voom_Mean_Variance_Trend.png",
  width = 2400,
  height = 1800,
  res = 300
)

v <- voom(
  dge_filtered,
  design,
  plot = TRUE
)

dev.off()

dim(v$E)

###Chunk 10: Fit the linear model###

fit <- lmFit(
  v,
  design
)

fit <- eBayes(
  fit
)

colnames(fit$coefficients)

###Chunk 11: Extract DE results###

results_all <- topTable(
  fit,
  coef = "ConditionUC",
  number = Inf,
  sort.by = "P"
)

head(results_all)
dim(results_all)

results_all$Gene <- rownames(results_all)

results_all <- results_all %>%
  select(Gene, everything())

write.csv(
  results_all,
  "output/DEG_All_Results_UC_vs_Healthy.csv",
  row.names = FALSE
)

###Chunk 12: Add DEG classification (Classify & summarize)###

results_all <- results_all %>%
  mutate(
    DEG_Status = case_when(
      adj.P.Val < 0.05 & logFC >= 1 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC <= -1 ~ "Downregulated",
      TRUE ~ "Not_Significant"
    )
  )

table(results_all$DEG_Status)

deg_summary <- results_all %>%
  count(DEG_Status)

deg_summary

results_all <- as.data.frame(results_all)

deg_summary <- dplyr::count(
  results_all,
  DEG_Status,
  name = "Number_of_Genes"
)

deg_summary

class(results_all)
str(results_all)

class(results_all$DEG_Status)
table(results_all$DEG_Status)

deg_summary <- as.data.frame(
  table(results_all$DEG_Status)
)

colnames(deg_summary) <- c(
  "DEG_Status",
  "Number_of_Genes"
)

deg_summary

write.csv(
  deg_summary,
  "output/DEG_Summary_UC_vs_Healthy.csv",
  row.names = FALSE
)

deg_summary

summary(results_all$adj.P.Val)

min(results_all$adj.P.Val, na.rm = TRUE)

summary(results_all$logFC)

range(results_all$logFC, na.rm = TRUE)

sum(results_all$adj.P.Val < 0.05, na.rm = TRUE)

sum(abs(results_all$logFC) >= 1, na.rm = TRUE)

sum(
  results_all$adj.P.Val < 0.05 &
    abs(results_all$logFC) >= 1,
  na.rm = TRUE
)

head(results_all[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], 10)

#25/08/2026

setwd("D:/BulkRNAseq/Deseq")

library(tidyverse)
library(readxl)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
library(ggrepel)

#SCRIPT A — Reload metadata and count data

metadata <- read_excel(
  "metadata.csv.xlsx",
  sheet = "data"
)

metadata <- as.data.frame(metadata)

merged <- read.delim(
  "ileum_raw_counts.csv",
  check.names = FALSE,
  row.names = 1
)

rownames(merged) <- make.unique(rownames(merged))

protein_coding_genes <- read_lines("pc_gene.txt")

merged <- merged[
  rownames(merged) %in% protein_coding_genes,
]

dim(merged)

#SCRIPT B — Rebuild metadata exactly as before

metadata <- metadata[
  metadata$Sample %in% colnames(merged),
]

metadata <- metadata[
  metadata$Condition %in% c("Healthy", "UC"),
]

metadata$Condition <- factor(metadata$Condition)
metadata$Condition <- droplevels(metadata$Condition)

metadata$Batch <- factor(metadata$Dataset)
metadata$Batch <- droplevels(metadata$Batch)

merged <- merged[
  ,
  colnames(merged) %in% metadata$Sample
]

metadata <- metadata[
  match(colnames(merged), metadata$Sample),
]

stopifnot(all(metadata$Sample == colnames(merged)))

dim(merged)
dim(metadata)

table(metadata$Batch, metadata$Condition)

#SCRIPT C — Clean any NA/invalid values

merged[is.na(merged)] <- 0
merged[!is.finite(as.matrix(merged))] <- 0

#SCRIPT D — Rebuild corrected log-CPM data

logcpm <- cpm(
  as.matrix(merged),
  log = TRUE,
  prior.count = 1
)

mod <- model.matrix(
  ~ Condition,
  data = metadata
)

corrected_data_limma_bc <- removeBatchEffect(
  logcpm,
  batch = metadata$Batch,
  design = mod
)

dim(corrected_data_limma_bc)
sum(is.na(corrected_data_limma_bc))

#SCRIPT E — Rebuild the DE analysis up to Script 10

dge <- DGEList(
  counts = merged,
  group = metadata$Condition
)

metadata$Condition <- relevel(
  factor(metadata$Condition),
  ref = "Healthy"
)

metadata$Batch <- factor(metadata$Batch)

design <- model.matrix(
  ~ Batch + Condition,
  data = metadata
)

colnames(design)

keep <- filterByExpr(
  dge,
  design = design
)

table(keep)

dge_filtered <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]
dim(dge_filtered)

#SCRIPT F — TMM normalization

dge_filtered <- calcNormFactors(
  dge_filtered,
  method = "TMM"
)

#SCRIPT G — Rebuild design and run voom

design <- model.matrix(
  ~ Batch + Condition,
  data = metadata
)

colnames(design)

qr(design)$rank
ncol(design)

v <- voom(
  dge_filtered,
  design,
  plot = TRUE
)

fit <- lmFit(
  v,
  design
)

fit <- eBayes(fit)

colnames(fit$coefficients)

results_all <- topTable(
  fit,
  coef = "ConditionUC",
  number = Inf,
  sort.by = "P"
)

results_all$Gene <- rownames(results_all)

results_all <- results_all %>%
  select(Gene, everything())

write.csv(
  results_all,
  "output/DEG_All_Results_UC_vs_Healthy.csv",
  row.names = FALSE
)

head(results_all)

#SCRIPT 11 — DEG classification

results_all <- results_all %>%
  mutate(
    DEG_Status = case_when(
      adj.P.Val < 0.05 & logFC >= 1 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC <= -1 ~ "Downregulated",
      TRUE ~ "Not_Significant"
    )
  )
table(results_all$DEG_Status)

deg_summary <- as.data.frame(
  table(results_all$DEG_Status)
)

colnames(deg_summary) <- c(
  "DEG_Status",
  "Number_of_Genes"
)

deg_summary

write.csv(
  deg_summary,
  "output/DEG_Summary_UC_vs_Healthy.csv",
  row.names = FALSE
)

#SCRIPT 12 — Check FDR results before proceeding

min(results_all$adj.P.Val, na.rm = TRUE)
sum(results_all$adj.P.Val < 0.05, na.rm = TRUE)
sum(results_all$P.Value < 0.05, na.rm = TRUE)
sum(
  results_all$P.Value < 0.05 &
    abs(results_all$logFC) >= 1,
  na.rm = TRUE
)
summary(results_all$adj.P.Val)

###Chunk 13: Significant & nominal DEG lists###

deg_significant <- results_all %>%
  filter(
    adj.P.Val < 0.05,
    abs(logFC) >= 1
  )
nrow(deg_significant)

write.csv(
  deg_significant,
  "output/Significant_DEGs_UC_vs_Healthy.csv",
  row.names = FALSE
)

###Chunk 14: MA plot & basic volcano plot###

png(
  "output/MA_Plot_UC_vs_Healthy.png",
  width = 2400,
  height = 1800,
  res = 300
)

plotMD(
  fit,
  column = "ConditionUC",
  main = "MA Plot: UC vs Healthy"
)

abline(
  h = c(-1, 1),
  lty = 2
)

dev.off()

#Volcano plot

results_volcano <- results_all %>%
  mutate(
    adj.P.Val = pmax(adj.P.Val, 1e-300),
    neg_log10_FDR = -log10(adj.P.Val)
  )

p_volcano <- ggplot(
  results_volcano,
  aes(
    x = logFC,
    y = neg_log10_FDR,
    color = DEG_Status
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: UC vs Healthy",
    x = "log2 Fold Change (UC / Healthy)",
    y = "-log10 Adjusted P-value",
    color = "DEG Status"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

print(p_volcano)

ggsave(
  "output/Volcano_Plot_UC_vs_Healthy.png",
  p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

min(results_all$adj.P.Val, na.rm = TRUE)

sum(results_all$adj.P.Val < 0.05, na.rm = TRUE)

sum(results_all$P.Value < 0.05, na.rm = TRUE)

sum(
  results_all$P.Value < 0.05 &
    abs(results_all$logFC) >= 1,
  na.rm = TRUE
)

nrow(deg_significant)

#Create an exploratory nominal DEG list

deg_nominal <- results_all %>%
  filter(
    P.Value < 0.05,
    abs(logFC) >= 1
  )
nrow(deg_nominal)

deg_nominal <- deg_nominal %>%
  mutate(
    Regulation = case_when(
      logFC >= 1 ~ "Upregulated_in_UC",
      logFC <= -1 ~ "Downregulated_in_UC"
    )
  )

table(deg_nominal$Regulation)

write.csv(
  deg_nominal,
  "output/Nominal_DEGs_Pvalue_0.05_LogFC_1_UC_vs_Healthy.csv",
  row.names = FALSE
)

#30/08/2026

setwd("D:/BulkRNAseq/Deseq")

library(tidyverse)
library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(factoextra)

dir.create("output", showWarnings = FALSE)

# Load metadata + raw counts, keep protein-coding genes

metadata <- read.csv("metadata.csv.xlsx")
merged   <- read.csv("ileum_raw_counts.csv", check.names = FALSE, row.names = 1)

rownames(merged) <- make.unique(rownames(merged))
protein_coding_genes <- read_lines("pc_gene.txt")
merged <- merged[rownames(merged) %in% protein_coding_genes, ]

metadata <- metadata %>% filter(Sample %in% colnames(merged))
metadata <- metadata[match(colnames(merged), metadata$Sample), ]
stopifnot(all(metadata$Sample == colnames(merged)))

str(metadata)

getwd()
list.files(pattern = "metadata")
readLines("metadata.csv.xlsx", n = 3)

rm(metadata)
metadata <- read.csv("metadata.csv.xlsx")
str(metadata)

readBin("metadata.csv", "raw", n = 4)

#30/08/2026 - SECOND TRIAL

setwd("D:/BulkRNAseq/Deseq")
getwd()
list.files()

if (!dir.exists("output")) dir.create("output")
library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(DESeq2)

library(readxl)

metadata <- read_excel("metadata.csv.xlsx", sheet = "data")
metadata <- as.data.frame(metadata)

merged <- read.delim("ileum_raw_counts.csv", check.names = FALSE, row.names = 1)

dim(metadata)
dim(merged)

rownames(merged) <- make.unique(rownames(merged))

protein_coding_file <- "pc_gene.txt"
protein_coding_genes <- read_lines(protein_coding_file)
merged <- merged[rownames(merged) %in% protein_coding_genes, ]

library(ggplot2)

plot_pca <- function(pca_df, var_exp, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"), color = "Condition") +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
}

plot_pcab <- function(pca_df, var_exp, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"), color = "Batch") +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
}

exists("pca_after_df")
exists("var_explained_after")

p1 <- plot_pcab(pca_after_df, var_explained_after, "PCA (log-CPM) After - by Batch")
ggsave("output/PCA_after_Batch.png", p1, width = 8, height = 6, dpi = 300)


setwd("D:/BulkRNAseq/Deseq")
list.files()

library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(ggplot2)

metadata <- read_excel("metadata.csv.xlsx", sheet = "data")   
metadata <- as.data.frame(metadata)

merged <- read.delim("ileum_raw_counts.csv", check.names = FALSE, row.names = 1)

rownames(merged) <- make.unique(rownames(merged))
protein_coding_genes <- read_lines("pc_gene.txt")
merged <- merged[rownames(merged) %in% protein_coding_genes, ]

metadata <- metadata[metadata$Sample %in% colnames(merged), ]
metadata <- metadata[match(colnames(merged), metadata$Sample), ]
stopifnot(all(metadata$Sample == colnames(merged)))

merged[is.na(merged)] <- 0
merged[!is.finite(as.matrix(merged))] <- 0

exclude <- c("GSE183620", "GSE193677")

metadata <- metadata[metadata$Condition %in% c("Healthy", "UC"), ]
metadata$Condition <- droplevels(as.factor(metadata$Condition))
metadata$Batch <- droplevels(as.factor(metadata$Dataset))

merged <- merged[, colnames(merged) %in% metadata$Sample]
metadata <- metadata[match(colnames(merged), metadata$Sample), ]
stopifnot(all(metadata$Sample == colnames(merged)))

exclude <- c("GSE183620", "GSE193677")

metadata <- metadata[metadata$Condition %in% c("Healthy", "UC"), ]
metadata$Condition <- droplevels(as.factor(metadata$Condition))
metadata$Batch <- droplevels(as.factor(metadata$Dataset))

merged <- merged[, colnames(merged) %in% metadata$Sample]
metadata <- metadata[match(colnames(merged), metadata$Sample), ]
stopifnot(all(metadata$Sample == colnames(merged)))

logcpm <- cpm(as.matrix(merged), log = TRUE, prior.count = 1)

mod <- model.matrix(~ Condition, data = metadata)
corrected_data_limma_bc <- removeBatchEffect(logcpm,
                                             batch = metadata$Batch,
                                             design = mod)
sum(is.na(corrected_data_limma_bc))

pca_before <- prcomp(t(logcpm), scale. = TRUE)
pca_before_df <- data.frame(pca_before$x, Condition = metadata$Condition, Batch = metadata$Batch)
var_explained_before <- round(100 * (pca_before$sdev^2 / sum(pca_before$sdev^2)), 2)

pca_after <- prcomp(t(corrected_data_limma_bc), scale. = TRUE)
pca_after_df <- data.frame(pca_after$x, Condition = metadata$Condition, Batch = metadata$Batch)
var_explained_after <- round(100 * (pca_after$sdev^2 / sum(pca_after$sdev^2)), 2)

plot_pca <- function(pca_df, var_exp, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"), color = "Condition") +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
}

plot_pcab <- function(pca_df, var_exp, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"), color = "Batch") +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
}

p1 <- plot_pcab(pca_after_df, var_explained_after, "PCA (log-CPM) After - by Batch")
ggsave("output/PCA_after_Batch.png", p1, width = 8, height = 6, dpi = 300)

var_explained_after

dim(merged)      # genes x samples going into logcpm
dim(metadata)   # should match ncol(merged)
nrow(logcpm)     # genes used for PCA
ncol(logcpm)     # samples used for PCA

length(protein_coding_genes)
nrow(merged)   # should be close to length(protein_coding_genes), not 46105 (the unfiltered count)

table(metadata$Batch, metadata$Condition)

sum(apply(logcpm, 1, var) == 0)
sum(apply(corrected_data_limma_bc, 1, var) == 0)

#differential expression analysis

library(tidyverse)
library(readxl)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
library(ggrepel)

write.csv(
  corrected_data_limma_bc,
  "output/logCPM_batch_corrected_limma.csv"
)

write.csv(
  metadata,
  "output/metadata_final_137_samples.csv",
  row.names = FALSE
)

write.csv(
  logcpm,
  "output/logCPM_before_batch_correction.csv"
)

lib_sizes <- colSums(merged)

summary(lib_sizes)

libsize_df <- data.frame(
  Sample = colnames(merged),
  Library_Size = lib_sizes,
  Condition = metadata$Condition,
  Batch = metadata$Batch
)

p_libsize <- ggplot(
  libsize_df,
  aes(x = Batch, y = Library_Size, fill = Condition)
) +
  geom_boxplot(alpha = 0.8) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Library Size Distribution by Batch",
    x = "Batch",
    y = "Total Raw Counts"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_libsize)

ggsave(
  "output/Library_Size_by_Batch.png",
  p_libsize,
  width = 9,
  height = 6,
  dpi = 300
)

dge <- DGEList(
  counts = merged,
  group = metadata$Condition
)
dge

#Filter lowly expressed genes

design <- model.matrix(
  ~ Batch + Condition,
  data = metadata
)

head(design)
dim(design)

keep <- filterByExpr(
  dge,
  design = design
)

table(keep)

dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

dim(dge_filtered)

#TMM normalization

dge_filtered <- calcNormFactors(
  dge_filtered,
  method = "TMM"
)

dge_filtered$samples


normalization_factors <- data.frame(
  Sample = rownames(dge_filtered$samples),
  dge_filtered$samples
)

write.csv(
  normalization_factors,
  "output/TMM_Normalization_Factors.csv",
  row.names = FALSE
)

#MDS plot after TMM normalization

png(
  "output/MDS_by_Condition.png",
  width = 2400,
  height = 1800,
  res = 300
)

plotMDS(
  dge_filtered,
  labels = metadata$Condition,
  col = ifelse(
    metadata$Condition == "Healthy",
    "tomato",
    "turquoise4"
  )
)

legend(
  "topright",
  legend = levels(metadata$Condition),
  col = c("tomato", "turquoise4"),
  pch = 16
)

dev.off()

#Recreate the design matrix carefully

metadata$Condition <- relevel(
  factor(metadata$Condition),
  ref = "Healthy"
)

metadata$Batch <- factor(metadata$Batch)

design <- model.matrix(
  ~ Batch + Condition,
  data = metadata
)

colnames(design)

qr(design)$rank
ncol(design)

#Run voom transformation

v <- voom(
  dge_filtered,
  design,
  plot = TRUE
)

png(
  "output/Voom_Mean_Variance_Trend.png",
  width = 2400,
  height = 1800,
  res = 300
)

v <- voom(
  dge_filtered,
  design,
  plot = TRUE
)

dev.off()

dim(v$E)

#Fit the linear model

fit <- lmFit(
  v,
  design
)

fit <- eBayes(
  fit
)

colnames(fit$coefficients)

#Extract DE results

results_all <- topTable(
  fit,
  coef = "ConditionUC",
  number = Inf,
  sort.by = "P"
)

head(results_all)
dim(results_all)

results_all$Gene <- rownames(results_all)

results_all <- results_all %>%
  select(Gene, everything())

write.csv(
  results_all,
  "output/DEG_All_Results_UC_vs_Healthy.csv",
  row.names = FALSE
)

#Add DEG classification

results_all <- results_all %>%
  mutate(
    DEG_Status = case_when(
      adj.P.Val < 0.05 & logFC >= 1 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC <= -1 ~ "Downregulated",
      TRUE ~ "Not_Significant"
    )
  )

table(results_all$DEG_Status)

deg_summary <- results_all %>%
  count(DEG_Status)

deg_summary

results_all <- as.data.frame(results_all)

deg_summary <- dplyr::count(
  results_all,
  DEG_Status,
  name = "Number_of_Genes"
)

deg_summary

class(results_all)
str(results_all)

class(results_all$DEG_Status)
table(results_all$DEG_Status)

deg_summary <- as.data.frame(
  table(results_all$DEG_Status)
)

colnames(deg_summary) <- c(
  "DEG_Status",
  "Number_of_Genes"
)

deg_summary

write.csv(
  deg_summary,
  "output/DEG_Summary_UC_vs_Healthy.csv",
  row.names = FALSE
)

deg_summary

summary(results_all$adj.P.Val)

min(results_all$adj.P.Val, na.rm = TRUE)

summary(results_all$logFC)

range(results_all$logFC, na.rm = TRUE)

sum(results_all$adj.P.Val < 0.05, na.rm = TRUE)

sum(abs(results_all$logFC) >= 1, na.rm = TRUE)

sum(
  results_all$adj.P.Val < 0.05 &
    abs(results_all$logFC) >= 1,
  na.rm = TRUE
)

head(results_all[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], 10)

#25/08/2026

setwd("D:/BulkRNAseq/Deseq")

library(tidyverse)
library(readxl)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
library(ggrepel)

#Reload metadata and count data

metadata <- read_excel(
  "metadata.csv.xlsx",
  sheet = "data"
)

metadata <- as.data.frame(metadata)

merged <- read.delim(
  "ileum_raw_counts.csv",
  check.names = FALSE,
  row.names = 1
)

rownames(merged) <- make.unique(rownames(merged))

protein_coding_genes <- read_lines("pc_gene.txt")

merged <- merged[
  rownames(merged) %in% protein_coding_genes,
]

dim(merged)

#Rebuild metadata exactly as before

metadata <- metadata[
  metadata$Sample %in% colnames(merged),
]

metadata <- metadata[
  metadata$Condition %in% c("Healthy", "UC"),
]

metadata$Condition <- factor(metadata$Condition)
metadata$Condition <- droplevels(metadata$Condition)

metadata$Batch <- factor(metadata$Dataset)
metadata$Batch <- droplevels(metadata$Batch)

merged <- merged[
  ,
  colnames(merged) %in% metadata$Sample
]

metadata <- metadata[
  match(colnames(merged), metadata$Sample),
]

stopifnot(all(metadata$Sample == colnames(merged)))

dim(merged)
dim(metadata)

table(metadata$Batch, metadata$Condition)

#Clean any NA/invalid values

merged[is.na(merged)] <- 0
merged[!is.finite(as.matrix(merged))] <- 0

#Rebuild corrected log-CPM data

logcpm <- cpm(
  as.matrix(merged),
  log = TRUE,
  prior.count = 1
)

mod <- model.matrix(
  ~ Condition,
  data = metadata
)

corrected_data_limma_bc <- removeBatchEffect(
  logcpm,
  batch = metadata$Batch,
  design = mod
)

dim(corrected_data_limma_bc)
sum(is.na(corrected_data_limma_bc))

#Rebuild the DE analysis up to Script 10

dge <- DGEList(
  counts = merged,
  group = metadata$Condition
)

metadata$Condition <- relevel(
  factor(metadata$Condition),
  ref = "Healthy"
)

metadata$Batch <- factor(metadata$Batch)

design <- model.matrix(
  ~ Batch + Condition,
  data = metadata
)

colnames(design)

keep <- filterByExpr(
  dge,
  design = design
)

table(keep)

dge_filtered <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]
dim(dge_filtered)

#TMM normalization

dge_filtered <- calcNormFactors(
  dge_filtered,
  method = "TMM"
)

#Rebuild design and run voom

design <- model.matrix(
  ~ Batch + Condition,
  data = metadata
)

colnames(design)

qr(design)$rank
ncol(design)

v <- voom(
  dge_filtered,
  design,
  plot = TRUE
)

fit <- lmFit(
  v,
  design
)

fit <- eBayes(fit)

colnames(fit$coefficients)

results_all <- topTable(
  fit,
  coef = "ConditionUC",
  number = Inf,
  sort.by = "P"
)

results_all$Gene <- rownames(results_all)

results_all <- results_all %>%
  select(Gene, everything())

write.csv(
  results_all,
  "output/DEG_All_Results_UC_vs_Healthy.csv",
  row.names = FALSE
)

head(results_all)

#SCRIPT 11 — DEG classification

results_all <- results_all %>%
  mutate(
    DEG_Status = case_when(
      adj.P.Val < 0.05 & logFC >= 1 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC <= -1 ~ "Downregulated",
      TRUE ~ "Not_Significant"
    )
  )
table(results_all$DEG_Status)

deg_summary <- as.data.frame(
  table(results_all$DEG_Status)
)

colnames(deg_summary) <- c(
  "DEG_Status",
  "Number_of_Genes"
)

deg_summary

write.csv(
  deg_summary,
  "output/DEG_Summary_UC_vs_Healthy.csv",
  row.names = FALSE
)

#Check FDR results before proceeding

min(results_all$adj.P.Val, na.rm = TRUE)
sum(results_all$adj.P.Val < 0.05, na.rm = TRUE)
sum(results_all$P.Value < 0.05, na.rm = TRUE)
sum(
  results_all$P.Value < 0.05 &
    abs(results_all$logFC) >= 1,
  na.rm = TRUE
)
summary(results_all$adj.P.Val)

#Extract significant DEGs

deg_significant <- results_all %>%
  filter(
    adj.P.Val < 0.05,
    abs(logFC) >= 1
  )
nrow(deg_significant)

write.csv(
  deg_significant,
  "output/Significant_DEGs_UC_vs_Healthy.csv",
  row.names = FALSE
)

#MA plot

png(
  "output/MA_Plot_UC_vs_Healthy.png",
  width = 2400,
  height = 1800,
  res = 300
)

plotMD(
  fit,
  column = "ConditionUC",
  main = "MA Plot: UC vs Healthy"
)

abline(
  h = c(-1, 1),
  lty = 2
)

dev.off()

#SCRIPT 15 — Volcano plot

results_volcano <- results_all %>%
  mutate(
    adj.P.Val = pmax(adj.P.Val, 1e-300),
    neg_log10_FDR = -log10(adj.P.Val)
  )

p_volcano <- ggplot(
  results_volcano,
  aes(
    x = logFC,
    y = neg_log10_FDR,
    color = DEG_Status
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: UC vs Healthy",
    x = "log2 Fold Change (UC / Healthy)",
    y = "-log10 Adjusted P-value",
    color = "DEG Status"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

print(p_volcano)

ggsave(
  "output/Volcano_Plot_UC_vs_Healthy.png",
  p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

min(results_all$adj.P.Val, na.rm = TRUE)

sum(results_all$adj.P.Val < 0.05, na.rm = TRUE)

sum(results_all$P.Value < 0.05, na.rm = TRUE)

sum(
  results_all$P.Value < 0.05 &
    abs(results_all$logFC) >= 1,
  na.rm = TRUE
)

nrow(deg_significant)

#Create an exploratory nominal DEG list

deg_nominal <- results_all %>%
  filter(
    P.Value < 0.05,
    abs(logFC) >= 1
  )
nrow(deg_nominal)

deg_nominal <- deg_nominal %>%
  mutate(
    Regulation = case_when(
      logFC >= 1 ~ "Upregulated_in_UC",
      logFC <= -1 ~ "Downregulated_in_UC"
    )
  )

table(deg_nominal$Regulation)

write.csv(
  deg_nominal,
  "output/Nominal_DEGs_Pvalue_0.05_LogFC_1_UC_vs_Healthy.csv",
  row.names = FALSE
)

#Contrast PCA (factoextra)

library(factoextra)

# Use voom-normalized expression matrix

pca_res <- prcomp(t(v$E), scale. = FALSE)

qq <- fviz_pca_ind(pca_res,
                   label = "none",
                   habillage = design$Condition,   # swap for your metadata's condition column
                   addEllipses = FALSE,
                   mean.point = FALSE,
                   title = "UC vs Healthy",
                   pointsize = 3) +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 11),
        axis.title = element_text(face = "plain", size = 10),
        axis.text = element_text(face = "plain", size = 9),
        legend.title = element_text(face = "plain", size = 10),
        legend.text = element_text(face = "plain", size = 9))

qq <- fviz_pca_ind(pca_res,
                   label = "none",
                   habillage = metadata$Condition,   # use metadata, not design
                   addEllipses = FALSE,
                   mean.point = FALSE,
                   title = "UC vs Healthy",
                   pointsize = 3) +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 11),
        axis.title = element_text(face = "plain", size = 10),
        axis.text = element_text(face = "plain", size = 9),
        legend.title = element_text(face = "plain", size = 10),
        legend.text = element_text(face = "plain", size = 9))

identical(colnames(v$E), metadata$Sample)

library(factoextra)


qq <- fviz_pca_ind(pca_res,
                   label = "none",
                   habillage = metadata$Condition,
                   addEllipses = FALSE,
                   mean.point = FALSE,
                   title = "UC vs Healthy",
                   pointsize = 3) +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 11),
        axis.title = element_text(face = "plain", size = 10),
        axis.text = element_text(face = "plain", size = 9),
        legend.title = element_text(face = "plain", size = 10),
        legend.text = element_text(face = "plain", size = 9))

png("output/ContrastPCA_UC_vs_Healthy.png", width = 7, height = 5, units = 'in', res = 300)
plot(qq)
dev.off()

###Chunk 15:  Refined volcano with highlighted genes###

library(ggrepel)

highlightgenes <- c("REG1A", "FDFT1", "SIK1", "DMBT1", "GGACT", "DUSP21")  # edit to your genes of interest

ls()

results_all <- results_all %>%
  mutate(
    Legend_Category = case_when(
      rownames(results_all) %in% highlightgenes & adj.P.Val < 0.05 ~ "Reported_in_Papers",
      adj.P.Val < 0.05 ~ "Significant",
      TRUE ~ "Not_significant"
    ),
    Legend_Category = factor(Legend_Category,
                             levels = c("Not_significant", "Significant", "Reported_in_Papers"))
  )

legend_colors <- c("Reported_in_Papers" = "red",
                   "Significant" = "#1e90ff",
                   "Not_significant" = "grey")

label_df <- results_all[rownames(results_all) %in% highlightgenes & results_all$adj.P.Val < 0.05, ]

volcano_refined <- ggplot(results_all, aes(x = logFC, y = -log10(adj.P.Val), color = Legend_Category)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = legend_colors) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = label_df, aes(label = rownames(label_df)),
                  color = "black", size = 2.5, box.padding = 0.4,
                  point.padding = 0.3, max.overlaps = 20, show.legend = FALSE) +
  theme_minimal() +
  labs(title = "UC vs Healthy", x = "log2 Fold Change",
       y = "-log10 Adjusted P value", color = "Legend") +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.title = element_text(size = 10), legend.text = element_text(size = 9),
        axis.title = element_text(size = 10), axis.text = element_text(size = 9))

ggsave("output/VolcanoPlot_Refined_UC_vs_Healthy.png", volcano_refined, width = 8, height = 6, dpi = 300)

###Chunk 16: Reactome GSEA###

library(clusterProfiler)

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
BiocManager::install("reactome.db")
BiocManager::install("ReactomePA")
library(ReactomePA)
BiocManager::install("enrichplot")
library(enrichplot)
BiocManager::install("stringr")
library(stringr)
BiocManager::install("forcats")
library(forcats)
library(clusterProfiler)

# Map gene symbols to ENTREZ IDs
ncbi_list <- bitr(rownames(results_all), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

library(dplyr)

results_all$names <- rownames(results_all)
deg_mapped <- results_all %>%
  left_join(ncbi_list, by = c("names" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

ngenes <- deg_mapped$logFC
names(ngenes) <- deg_mapped$ENTREZID
ngenes <- sort(ngenes, decreasing = TRUE)

enp_gsea <- gsePathway(geneList = ngenes, organism = "human", verbose = FALSE)

# Ranked gene list by logFC

ngenes <- deg_mapped$logFC
names(ngenes) <- deg_mapped$ENTREZID
ngenes <- sort(ngenes, decreasing = TRUE)

enp_gsea <- gsePathway(geneList = ngenes, organism = "human", verbose = FALSE)
enp_gsea <- setReadable(enp_gsea, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

pathways <- enp_gsea@result %>% arrange(p.adjust)
top_pathways <- pathways %>% arrange(desc(abs(NES)))
top_pathways$Description <- str_wrap(top_pathways$Description, width = 50)
top_pathways <- top_pathways[1:20, ] %>% mutate(Description = fct_reorder(Description, NES))

library(ggplot2)

pathway_plot <- ggplot(top_pathways, aes(x = NES, y = Description, color = p.adjust, size = setSize)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#0072B2", high = "#D55E00", name = "FDR (p.adjust)") +
  scale_size(range = c(3, 10), name = "Gene Set Size") +
  labs(title = "Top 20 Reactome Pathways", x = "Normalized Enrichment Score", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 15, hjust = 0.5),
        legend.title = element_text(size = 12), legend.text = element_text(size = 11))

ggsave("output/PathwaysPlot_UC_vs_Healthy.png", pathway_plot, width = 8, height = 6)
write.csv(pathways, "output/Reactome_Pathways_UC_vs_Healthy.csv")

## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)

library(hdf5r)

library(multtest)
library(metap)

library(SingleR)
library('scuttle')

Ov.integrated = readRDS('Ov.integrated.rds')
# Extract raw counts and metadata to create SingleCellExperiment object
counts <- Ov.integrated@assays$RNA@counts 
counts

metadata <- Ov.integrated@meta.data
metadata

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(Ov.integrated@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

## create sample id
colData(sce)$sample_id = paste0(colData(sce)$age, '_', colData(sce)$hash.ID)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "sample_id")]
# groups$cluster_id = gsub('_', '', groups$cluster_id)
# groups$sample_id = gsub('_', '', groups$sample_id)
groups$cluster_id = gsub('_', '\\.', groups$cluster_id)
groups$sample_id = gsub('_', '', groups$sample_id)
groups

dim(counts((sce)))
(counts((sce)))
colData((sce))

# Aggregate across cluster-sample groups
library(Matrix.utils)
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

dim(pb)
pb[1:10, 1:6]

# splitf = sapply(strsplit(rownames(pb), '_'), function(x)x[1]) 
splitf = sapply(strsplit(rownames(pb), '\\.'), function(x)x[1]) 

split.data.frame(pb, factor(splitf))

library(magrittr)
pb = split.data.frame(pb, 
                 factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                  # stringr::str_extract(rownames(pb), "months.*?HTO.*?$")))
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
                  # stringr::str_extract(rownames(pb), "HTO.")))
                  # rownames(pb)))

str(pb)

# Print out the table of cells in each cluster-sample group
table(sce$cluster_id, sce$sample_id)

pb$B11%>%head()

metadata = colData(sce)[, c('celltype', 'sample_id', 'age')]
metadata$sample_id = gsub('_', '', metadata$sample_id)
metadata
metadata = metadata%>%unique()
metadata
colnames(metadata) = c('cluster_id', 'sample_id', 'group_id')
metadata$group_id

clusters <- levels(metadata$cluster_id)
clusters

# Subset the metadata to only the Mac0 cells
cluster_ind = 1
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[cluster_ind]), ]
cluster_metadata$sample_id = cluster_metadata$sample_id%>%as.factor()
cluster_metadata$group_id = cluster_metadata$group_id%>%as.factor()
head(cluster_metadata)

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id
head(cluster_metadata)

# Subset the counts to only the B cells
counts <- pb[[clusters[1]]]
clusters
colnames(counts)
rownames(cluster_metadata)

counts%>%head()
counts

# cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
counts
cluster_counts = data.frame(counts[, match( rownames(cluster_metadata), colnames(counts))])
cluster_counts%>%head()

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))     

### DESeq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ group_id)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
DESeq2::plotPCA(rld, intgroup = "group_id")

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
rld_cor

# Plot heatmap
library(pheatmap)
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(nrow(rld_cor)/2), colorRampPalette(colors = c( "orange", "red", "purple"))(nrow(rld_cor)/2))
pheatmap(rld_cor, annotation = data.frame(Var =(cluster_metadata[,'group_id', drop = F])), color = my.colors)

as.factor(as.vector(cluster_metadata$group_id))

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Output results of Wald test for contrast for stim vs ctrl
levels(cluster_metadata$group_id)[2]
cluster_metadata$group_id
levels(cluster_metadata$group_id)[1]

contrast <- c("group_id", levels(cluster_metadata$group_id)[2], levels(cluster_metadata$group_id)[1])

# resultsNames(dds)
res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

res <- lfcShrink(dds, 
                 contrast =  contrast,
                 type = 'ashr',
                 res=res)

res[order(res$log2FoldChange, decreasing = T),]
res[grep('Cxcl9', rownames(res)),]
# write.csv(res[order(res$log2FoldChange, decreasing = T),], 'Mac0_pseudoBulkDE.csv')

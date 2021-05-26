## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
# BiocManager::install("hdf5r")
library(hdf5r)
library(dplyr)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
FeatureBC3 <- Read10X_h5("./ov3bc/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
Ov3.hto <- FeatureBC3$`Antibody Capture`
Ov3.hto <- as.data.frame(Ov3.hto)
Ov3.umi <- FeatureBC3$`Gene Expression`
HTOs <- c("HTO0301","HTO0302","HTO0303")
rownames(Ov3.hto) <- HTOs

# Select cell barcodes detected by both RNA and HTO. 
joint.bcs <- intersect(colnames(Ov3.umi), colnames(Ov3.hto))
Ov3.umi <- Ov3.umi[,joint.bcs]
Ov3.hto <- Ov3.hto[,joint.bcs]




## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov3.hashtag <- CreateSeuratObject(counts = Ov3.umi)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov3.hashtag[["HTO"]] <- CreateAssayObject(counts = Ov3.hto)


#QC removing high mt count cells
Ov3.hashtag[["percent.mt"]] <- PercentageFeatureSet(Ov3.hashtag, pattern = "^mt-")
head(Ov3.hashtag@meta.data, 5)

VlnPlot(Ov3.hashtag, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Ov3.hashtag <- subset(Ov3.hashtag, subset = nFeature_RNA > 300 & nFeature_RNA < 5500 & percent.mt < 5)
Ov3.hashtag

#Data normalization with log normalization
Ov3.hashtag <- NormalizeData(Ov3.hashtag)
Ov3.hashtag <- NormalizeData(Ov3.hashtag, assay = "HTO", normalization.method = "CLR")
#Find and scale variable features
Ov3.hashtag <- FindVariableFeatures(Ov3.hashtag, selection.method = "mean.var.plot")
Ov3.hashtag <- ScaleData(Ov3.hashtag, features = VariableFeatures(Ov3.hashtag))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov3.hashtag <- HTODemux(Ov3.hashtag, assay = "HTO", positive.quantile = 0.99)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
table(Ov3.hashtag$HTO_classification.global)

Idents(Ov3.hashtag) <- "HTO_maxID"
RidgePlot(Ov3.hashtag, assay = "HTO", features = rownames(Ov3.hashtag[["HTO"]]), ncol = 3)



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
FeatureScatter(Ov3.hashtag, feature1 = "hto_HTO0301", feature2 = "hto_HTO0302")
FeatureScatter(Ov3.hashtag, feature1 = "hto_HTO0301", feature2 = "hto_HTO0303")
FeatureScatter(Ov3.hashtag, feature1 = "hto_HTO0303", feature2 = "hto_HTO0302")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Idents(Ov3.hashtag) <- "HTO_classification.global"
VlnPlot(Ov3.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# # Remove the negative cells
# Ov3.hashtag.subset <- subset(Ov3.hashtag, idents = "Negative", invert = TRUE)
# 
# # Run PCA on most variable features
# Ov3.hashtag <- FindVariableFeatures(Ov3.hashtag, selection.method = "mean.var.plot")
# Ov3.hashtag <- ScaleData(Ov3.hashtag, features = VariableFeatures(Ov3.hashtag))
# Ov3.hashtag <- RunPCA(Ov3.hashtag)
# Ov3.hashtag <- RunTSNE(Ov3.hashtag, dims = 1:5, perplexity = 100)
# DimPlot(Ov3.hashtag) + NoLegend()


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#Removing negative cells
Ov3.hashtag.subset <- subset(Ov3.hashtag, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(Ov3.hashtag.subset) <- "HTO"
Ov3.hashtag.subset <- ScaleData(Ov3.hashtag.subset, features = rownames(Ov3.hashtag.subset), 
                                 verbose = FALSE)
Ov3.hashtag.subset <- RunPCA(Ov3.hashtag.subset, features = rownames(Ov3.hashtag.subset), approx = FALSE)
Ov3.hashtag.subset <- RunTSNE(Ov3.hashtag.subset, dims = 1:3, perplexity = 100, check_duplicates = FALSE)
DimPlot(Ov3.hashtag.subset)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
HTOHeatmap(Ov3.hashtag, assay = "HTO", ncells = 5000)

# Extract the singlets
Ov3.singlet <- subset(Ov3.hashtag, idents = "Singlet")

# Select the top 1000 most variable features
Ov3.singlet <- FindVariableFeatures(Ov3.singlet, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
Ov3.singlet <- ScaleData(Ov3.singlet, features = VariableFeatures(Ov3.singlet))

# Run PCA
Ov3.singlet <- RunPCA(Ov3.singlet, features = VariableFeatures(Ov3.singlet))

# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
Ov3.singlet <- FindNeighbors(Ov3.singlet, reduction = "pca", dims = 1:10)
Ov3.singlet <- FindClusters(Ov3.singlet, resolution = 0.6, verbose = FALSE)
Ov3.singlet <- RunTSNE(Ov3.singlet, reduction = "pca", dims = 1:10)

# Projecting singlet identities on TSNE visualization
DimPlot(Ov3.singlet, group.by = "HTO_classification")



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# HTOHeatmap(Ov3.hashtag, assay = "HTO", ncells = 3000)
# ```
# Cluster visualizing cells with usual scRNAseq workflow
# 
# ```{r}
# #Extracting singlets
# Ov3.singlet <- subset(Ov3.hashtag, idents = "Singlet")
# 
# 
# #Selecting top 1000 variable features
# Ov3.singlet <- FindVariableFeatures(Ov3.singlet, selection.method = "mean.var.plot")
# 
# #Scaling RNA data
# Ov3.singlet <- ScaleData(Ov3.singlet, features = VariableFeatures(Ov3.singlet))
# 
# #Run PCA
# Ov3.singlet <- RunPCA(Ov3.singlet, features = VariableFeatures(Ov3.singlet))
# 




## -------------------------------------------------------------------------------------------------------------------------------------------------------------
DimPlot(Ov3.singlet)

DimHeatmap(Ov3.singlet)
DimHeatmap(Ov3.singlet, dims = 1:15, cells = 1000, balanced  = T)



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
print(Ov3.singlet[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Ov3.singlet, dim = 1:2, reduction = "pca")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov3.singlet <- JackStraw(Ov3.singlet, num.replicate = 100)
Ov3.singlet <- ScoreJackStraw(Ov3.singlet, dims = 1:20)
JackStrawPlot(Ov3.singlet, dims = 1:20)
ElbowPlot(Ov3.singlet)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov3.singlet <- FindNeighbors(Ov3.singlet, dims = 1:13)
Ov3.singlet <- FindClusters(Ov3.singlet, resolution = 0.5)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov3.singlet <- RunUMAP(Ov3.singlet, dims = 1:13)

pdf(file = "ov3HTO_UMAP_clustering_v1.pdf")
DimPlot(Ov3.singlet, reduction = "umap")
dev.off()


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(Ov3.singlet, file = "./Ov3HTOumap.rds")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#identify  DE markers for all lusters compare to remaining of the cells all at once
#v4: return.thresh = 0.05, filtering out genes with p_val_adj >0.05
Ov3singlet.markers <- FindAllMarkers(Ov3.singlet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, return.thresh = 0.05)
Ov3singlet.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

write.csv(Ov3singlet.markers, file = "ov3HTOmarkersv1.csv")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
 pdf(file = "Ov3HTO_cluster_DEFeature_Heatmap_v1.pdf", height = 10)
top10 <- Ov3singlet.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Ov3.singlet,features = top10$gene) + NoLegend()
dev.off()



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
library(SingleR)
library('scuttle')
Ov3singlet.counts <- GetAssayData(Ov3.singlet)
singler3.singlet <- SingleR(test = Ov3singlet.counts, ref = ref.se, labels = ref.se$label.main, clusters = Ov3.singlet$seurat_clusters )
singler3.singlet



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
table(singler3.singlet$labels)
Ov3HTO_v1_annotation <- data.frame(singler3.singlet)

write.csv(Ov3HTO_v1_annotation, file = "ov3HTO_cluster_annotation.csv")



## -------------------------------------------------------------------------------------------------------------------------------------------------------------

pdf(file = "OV3singlet_UMAP_HTO_distribution.pdf")
DimPlot(Ov3.singlet, group.by = "HTO_classification")
FeaturePlot(Ov3.singlet, features = "HTO0301")
FeaturePlot(Ov3.singlet, features = "HTO0302")
FeaturePlot(Ov3.singlet, features = "HTO0303")
dev.off()


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
new.cluster.ids <- c("Mac1", "Mac2", "Mac3", "NK", "Mono1", "Mono2", "Mac4", "Neutrophil", "DC1", "Fibroblast", "Mac5", "DC2","endothelial")
 names(new.cluster.ids) <- levels(Ov3.singlet)
Ov3.singlet <- RenameIdents(Ov3.singlet, new.cluster.ids)
pdf(file = "Ov3singlet_UMAP_v1_annotated.pdf")
DimPlot(Ov3.singlet, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
FeaturePlot(Ov3.singlet, features = "Mrc1")
FeaturePlot(Ov3.singlet, features = "Msr1")
FeaturePlot(Ov3.singlet, features = "Arg1")
FeaturePlot(Ov3.singlet, features = "C1qa")
FeaturePlot(Ov3.singlet, features = "C1qc")

FeaturePlot(Ov3.singlet, features = "Itgax")
FeaturePlot(Ov3.singlet, features = "Cd86")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
FeatureBC15 <- Read10X_h5("./ov15bc/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
Ov15.hto <- FeatureBC15$`Antibody Capture`
Ov15.hto <- as.data.frame(Ov15.hto)
Ov15.umi <- FeatureBC15$`Gene Expression`
HTOs <- c("HTO0302","HTO0303","HTO0301")
rownames(Ov15.hto) <- HTOs

# Select cell barcodes detected by both RNA and HTO. 
Ov15joint.bcs <- intersect(colnames(Ov15.umi), colnames(Ov15.hto))
Ov15.umi <- Ov15.umi[,Ov15joint.bcs]
Ov15.hto <- Ov15.hto[,Ov15joint.bcs]




## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov15.hashtag <- CreateSeuratObject(counts = Ov15.umi)





## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov15.hashtag[["HTO"]] <- CreateAssayObject(counts = Ov15.hto)


#QC removing high mt count cells
Ov15.hashtag[["percent.mt"]] <- PercentageFeatureSet(Ov15.hashtag, pattern = "^mt-")
head(Ov15.hashtag@meta.data, 5)

VlnPlot(Ov15.hashtag, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Ov15.hashtag <- subset(Ov15.hashtag, subset = nFeature_RNA > 300 & nFeature_RNA < 5500 & percent.mt < 5)
Ov15.hashtag

#Data normalization with log normalization
Ov15.hashtag <- NormalizeData(Ov15.hashtag)
Ov15.hashtag <- NormalizeData(Ov15.hashtag, assay = "HTO", normalization.method = "CLR")
#Find and scale variable features
Ov15.hashtag <- FindVariableFeatures(Ov15.hashtag, selection.method = "mean.var.plot")
Ov15.hashtag <- ScaleData(Ov15.hashtag, features = VariableFeatures(Ov15.hashtag))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov15.hashtag <- HTODemux(Ov15.hashtag, assay = "HTO", positive.quantile = 0.99)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
table(Ov15.hashtag$HTO_classification.global)

Idents(Ov15.hashtag) <- "HTO_maxID"
RidgePlot(Ov15.hashtag, assay = "HTO", features = rownames(Ov15.hashtag[["HTO"]]), ncol = 3)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
FeatureScatter(Ov15.hashtag, feature1 = "hto_HTO0301", feature2 = "hto_HTO0302")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Idents(Ov15.hashtag) <- "HTO_classification.global"
VlnPlot(Ov15.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#Removing negative cells
Ov15.hashtag.subset <- subset(Ov15.hashtag, idents = "Negative", invert = TRUE)

#Calculating distant matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = Ov15.hashtag.subset, assay = "HTO"))))

#Calculating tSNE embedding with a distance matrix
Ov15.hashtag.subset <- RunTSNE(Ov15.hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(Ov15.hashtag.subset)



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
HTOHeatmap(Ov15.hashtag, assay = "HTO", ncells = 3000)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#Extracting singlets
Ov15.singlet <- subset(Ov15.hashtag, idents = "Singlet")


#Selecting top 1000 variable features
Ov15.singlet <- FindVariableFeatures(Ov15.singlet, selection.method = "mean.var.plot")

#Scaling RNA data
Ov15.singlet <- ScaleData(Ov15.singlet, features = VariableFeatures(Ov15.singlet))

#Run PCA
Ov15.singlet <- RunPCA(Ov15.singlet, features = VariableFeatures(Ov15.singlet))





## -------------------------------------------------------------------------------------------------------------------------------------------------------------
DimPlot(Ov15.singlet)

DimHeatmap(Ov15.singlet) 
DimHeatmap(Ov15.singlet, dims = 1:15, cells = 1000, balanced  = T)



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
print(Ov15.singlet[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Ov15.singlet, dim = 1:2, reduction = "pca")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov15.singlet <- JackStraw(Ov15.singlet, num.replicate = 100)
Ov15.singlet <- ScoreJackStraw(Ov15.singlet, dims = 1:20)
JackStrawPlot(Ov15.singlet, dims = 1:20)
ElbowPlot(Ov15.singlet)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov15.singlet <- FindNeighbors(Ov15.singlet, dims = 1:14)
Ov15.singlet <- FindClusters(Ov15.singlet, resolution = 0.5)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov15.singlet <- RunUMAP(Ov15.singlet, dims = 1:14)

pdf(file = "ov15HTO_UMAP_clustering_v1.pdf")
DimPlot(Ov15.singlet, reduction = "umap")
dev.off()


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(Ov15.singlet, file = "./Ov15HTOumap.rds")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
#identify  DE markers for all lusters compare to remaining of the cells all at once
#v4: return.thresh = 0.05, filtering out genes with p_val_adj >0.05
Ov15singlet.markers <- FindAllMarkers(Ov15.singlet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, return.thresh = 0.05)
Ov15singlet.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

write.csv(Ov15singlet.markers, file = "ov15HTOmarkersv1.csv")



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
 pdf(file = "Ov15HTO_cluster_DEFeature_Heatmap_v1.pdf", height = 10)
top10 <- Ov15singlet.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Ov15.singlet,features = top10$gene) + NoLegend()
dev.off()



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
library(SingleR)
library('scuttle')
Ov15singlet.counts <- GetAssayData(Ov15.singlet)
singler15.singlet <- SingleR(test = Ov15singlet.counts, ref = ref.se, labels = ref.se$label.main, clusters = Ov15.singlet$seurat_clusters )
singler15.singlet



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
table(singler15.singlet$labels)
Ov15HTO_v1_annotation <- data.frame(singler15.singlet)

write.csv(Ov15HTO_v1_annotation, file = "ov15HTO_cluster_annotation.csv")



## -------------------------------------------------------------------------------------------------------------------------------------------------------------

pdf(file = "OV15singlet_UMAP_HTO_distribution.pdf")
DimPlot(Ov15.singlet, group.by = "HTO_classification")
FeaturePlot(Ov15.singlet, features = "HTO0301")
FeaturePlot(Ov15.singlet, features = "HTO0302")
FeaturePlot(Ov15.singlet, features = "HTO0303")
dev.off()


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
new.cluster.ids <- c("Mac1", "DC", "Mac2", "Mono1", "Neutrophil", "Mono2", "NK", "B", "Mac3", "fibroblast", "Mac4", "CD8+ T cells","endothelial", "Mac5")
 names(new.cluster.ids) <- levels(Ov15.singlet)
Ov15.singlet <- RenameIdents(Ov15.singlet, new.cluster.ids)
pdf(file = "Ov15singlet_UMAP_v1_annotated.pdf")
DimPlot(Ov15.singlet, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

## -------------------------------------------------------------------------------------------------------------------------------------------------------------
FeaturePlot(Ov15.singlet, features = "Mrc1")
FeaturePlot(Ov15.singlet, features = "Msr1")
FeaturePlot(Ov15.singlet, features = "Cd163")
FeaturePlot(Ov15.singlet, features = "Arg1")
FeaturePlot(Ov15.singlet, features = "C1qa")
FeaturePlot(Ov15.singlet, features = "C1qb")
FeaturePlot(Ov15.singlet, features = "C1qc")
FeaturePlot(Ov15.singlet, features = "Fam20c")
FeaturePlot(Ov15.singlet, features = "Cstb")
FeaturePlot(Ov15.singlet, features = "Cd86")
FeaturePlot(Ov15.singlet, features = "Itgax")


## ----pressure, echo=FALSE-------------------------------------------------------------------------------------------------------------------------------------
plot(pressure)


## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)

library(hdf5r)

#BiocManager::install('multtest') 
#install.packages('metap')

library(multtest)
library(metap)

library(SingleR)
library('scuttle')



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
FeatureBC3 <- Read10X_h5("myeloid_data/ov3bc/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
Ov3.hto <- FeatureBC3$`Antibody Capture`
Ov3.hto <- as.data.frame(Ov3.hto)
Ov3.umi <- FeatureBC3$`Gene Expression`
HTOs3 <- c("HTO0301","HTO0302","HTO0303")
rownames(Ov3.hto) <- HTOs3

# Select cell barcodes detected by both RNA and HTO. 
joint.bcs3 <- intersect(colnames(Ov3.umi), colnames(Ov3.hto))
Ov3.umi <- Ov3.umi[,joint.bcs3]
Ov3.hto <- Ov3.hto[,joint.bcs3]

Ov3.umi


FeatureBC15 <- Read10X_h5("myeloid_data/ov15bc/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
Ov15.hto <- FeatureBC15$`Antibody Capture`
Ov15.hto <- as.data.frame(Ov15.hto)
Ov15.umi <- FeatureBC15$`Gene Expression`
HTOs15 <- c("HTO0302","HTO0303","HTO0301")
rownames(Ov15.hto) <- HTOs15

# Select cell barcodes detected by both RNA and HTO. 
Ov15joint.bcs <- intersect(colnames(Ov15.umi), colnames(Ov15.hto))
Ov15.umi <- Ov15.umi[,Ov15joint.bcs]
Ov15.hto <- Ov15.hto[,Ov15joint.bcs]


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov3.hashtag <- CreateSeuratObject(counts = Ov3.umi, project ="ovmonths3")

Ov3.hashtag[["HTO"]] <- CreateAssayObject(counts = Ov3.hto)
Ov3.hashtag[["age"]] <- "months3"

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Ov3.hashtag <- NormalizeData(Ov3.hashtag, assay = "HTO", normalization.method = "CLR")
GetAssayData(Ov3.hashtag)

#demultiplex HTO
Ov3.hashtag <- HTODemux(Ov3.hashtag, assay = "HTO", positive.quantile = 0.99)
GetAssayData(Ov3.hashtag)
# Global classification results
table(Ov3.hashtag$HTO_classification.global)
Idents(Ov3.hashtag) <- "HTO_classification.global"
#Extracting singlets
Ov3.singlet <- subset(Ov3.hashtag, idents = "Singlet")
GetAssayData(Ov3.singlet)


Ov15.hashtag <- CreateSeuratObject(counts = Ov15.umi, project ="ovmonths15")
Ov15.hashtag[["HTO"]] <- CreateAssayObject(counts = Ov15.hto)
Ov15.hashtag[["age"]] <- "months15"

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Ov15.hashtag <- NormalizeData(Ov15.hashtag, assay = "HTO", normalization.method = "CLR")

#demultiplex HTO
Ov15.hashtag <- HTODemux(Ov15.hashtag, assay = "HTO", positive.quantile = 0.99)
# Global classification results
table(Ov15.hashtag$HTO_classification.global)
Idents(Ov15.hashtag) <- "HTO_classification.global"
#Extracting singlets
Ov15.singlet <- subset(Ov15.hashtag, idents = "Singlet")



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov.combined <- merge(Ov3.singlet, y = Ov15.singlet, add.cell.ids = c("months3", "months15"))
GetAssayData(Ov.combined)
GetAssayData(Ov3.singlet)




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#QC removing high mt count cells
Ov.combined[["percent.mt"]] <- PercentageFeatureSet(Ov.combined, pattern = "^mt-")
head(Ov.combined@meta.data, 5)

VlnPlot(Ov.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Ov.combined <- subset(Ov.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 5)
GetAssayData(Ov.combined)





## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov.list <- SplitObject(Ov.combined, split.by = "age")
GetAssayData(Ov.list$months15)
Ov.list <- lapply(X = Ov.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

Ov.anchors <- FindIntegrationAnchors(object.list = Ov.list, dims = 1:20)
Ov.integrated <- IntegrateData(anchorset = Ov.anchors, dims = 1:20)
GetAssayData(Ov.integrated)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DefaultAssay(Ov.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
Ov.integrated <- ScaleData(Ov.integrated, verbose = FALSE)
Ov.integrated <- RunPCA(Ov.integrated, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
Ov.integrated <- RunUMAP(Ov.integrated, reduction = "pca", dims = 1:20)
Ov.integrated <- FindNeighbors(Ov.integrated, reduction = "pca", dims = 1:20)
Ov.integrated <- FindClusters(Ov.integrated, resolution = 0.5)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p1 <- DimPlot(Ov.integrated, reduction = "umap", group.by = "age")
p2 <- DimPlot(Ov.integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
p1
p2


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pdf("OvIntegrated_UMAP_splitbyage.pdf", width = 15)
DimPlot(Ov.integrated, reduction = "umap", split.by = "age")
dev.off()


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#identify  DE markers for all lusters compare to remaining of the cells all at once
#v4: return.thresh = 0.05, filtering out genes with p_val_adj >0.05
# Ov.integrated.markers <- FindAllMarkers(Ov.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, return.thresh = 0.05)
# Ov.integrated.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# 
# write.csv(Ov.integrated.markers, file = "OvIntegrated_markersv1.csv")



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DefaultAssay(Ov.integrated) <- "RNA"
# mac0.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 0, grouping.var = "age", verbose = FALSE)
# head(mac0.markers)
# mono1.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 1, grouping.var = "age", verbose = FALSE)
# mac2.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 2, grouping.var = "age", verbose = FALSE)
# neurophil.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 3, grouping.var = "age", verbose = FALSE)
# mono4.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 4, grouping.var = "age", verbose = FALSE)
# head(mono4.markers)
# DC5.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 5, grouping.var = "age", verbose = FALSE)
# NK6.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 6, grouping.var = "age", verbose = FALSE)
# mac7.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 7, grouping.var = "age", verbose = FALSE)
# T8.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 8, grouping.var = "age", verbose = FALSE)
# DC9.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 9, grouping.var = "age", verbose = FALSE)
# fibroblast10.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 10, grouping.var = "age", verbose = FALSE)
# B11.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 11, grouping.var = "age", verbose = FALSE)
# Mac12.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 12, grouping.var = "age", verbose = FALSE)
# Mac13.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 13, grouping.var = "age", verbose = FALSE)
# DC14.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 14, grouping.var = "age", verbose = FALSE)
# endothelial15.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 15, grouping.var = "age", verbose = FALSE)
# Mac16.markers <- FindConservedMarkers(Ov.integrated, ident.1 = 16, grouping.var = "age", verbose = FALSE)





## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov.integrated <- RenameIdents(Ov.integrated, `0` = "Mac0", `1` = "Mono1", `2` = "Mac2", 
    `3` = "Neutrophil", `4` = "Mono4", `5` = "DC5", `6` = "NK6", `7` = "Mac7", `8` = "T8", `9` = "DC9", 
    `10` = "Fibroblast10", `11` = "B11", `12` = "Mac12", `13` = "Mac13", `14` = "DC14", `15` = "Endothelial15", `16` = "Mac16")
pdf(file = "./OvIntegrated_UMAP_annotated.pdf")
DimPlot(Ov.integrated, label = TRUE)
dev.off()


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Idents(Ov.integrated) <- factor(Idents(Ov.integrated), levels = c("Mac0", "Mono1","Mac2","Neutrophil", "Mono4","DC5",  "NK6",  "Mac7", "T8", "DC9", "Fibroblast10","B11", "Mac12",  "Mac13", "DC14",  "Endothelial15", "Mac16"))

markers.to.plot <- c("C1qa", "Cd68", "Mrc1",  "Ly6c2", "Ccr2", "Arg1", "Spp1",  "S100a9", "S100a8","Treml1", "Ly6i", "Cd209a", "Clec10a", "Nkg7", "Klre1", "Mki67", "Cd300e", "Xcl1", "Cd3d", "Xcr1", "Clec9a",  "Col4a1", "Col1a2", "Ly6d", "Ms4a1", "Asb2", "Itgax",  
    "Il4i1", "Ccr7","Cldn5", "Cdh5", "Ccl24")
pdf(file ="./OvIntegrated_Markerexpression.pdf", height = 10, width = 15)
DotPlot(Ov.integrated, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
    split.by = "age") + RotatedAxis()
dev.off()


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Ov.integrated$celltype.age <- paste(Idents(Ov.integrated), Ov.integrated$age, sep = "_")
Ov.integrated$celltype <- Idents(Ov.integrated)
Idents(Ov.integrated) <- "celltype.age"
Mac0.aging <- FindMarkers(Ov.integrated, ident.1 = "Mac0_months3", ident.2 = "Mac0_months15", verbose = FALSE)
head(Mac0.aging, n = 15)
write.csv(Mac0.aging, "Mac0_aging.csv")

Mono1.aging <- FindMarkers(Ov.integrated, ident.1 = "Mono1_months3", ident.2 = "Mono1_months15", verbose = FALSE)
head(Mono1.aging, n = 15)
write.csv(Mono1.aging, "Mono1_aging.csv")

Mac2.aging <- FindMarkers(Ov.integrated, ident.1 = "Mac2_months3", ident.2 = "Mac2_months15", verbose = FALSE)
head(Mac2.aging, n = 15)
write.csv(Mac2.aging, "Mac2_aging.csv")

Neutrophil.aging <- FindMarkers(Ov.integrated, ident.1 = "Neutrophil_months3", ident.2 = "Neutrophil_months15", verbose = FALSE)
head(Neutrophil.aging, n = 15)
write.csv(Neutrophil.aging, "Neutrophil_aging.csv")

Mono4.aging <- FindMarkers(Ov.integrated, ident.1 = "Mono4_months3", ident.2 = "Mono4_months15", verbose = FALSE)
head(Mono4.aging, n = 15)
write.csv(Mono4.aging, "Mono4_aging.csv")

DC5.aging <- FindMarkers(Ov.integrated, ident.1 = "DC5_months3", ident.2 = "DC5_months15", verbose = FALSE)
head(DC5.aging, n = 15)
write.csv(DC5.aging, "DC5_aging.csv")

DC9.aging <- FindMarkers(Ov.integrated, ident.1 = "DC9_months3", ident.2 = "DC9_months15", verbose = FALSE)
head(DC9.aging, n = 15)
write.csv(DC9.aging, "DC9_aging.csv")

DC14.aging <- FindMarkers(Ov.integrated, ident.1 = "DC14_months3", ident.2 = "DC14_months15", verbose = FALSE)
head(DC14.aging, n = 15)
write.csv(DC14.aging, "DC14_aging.csv")

Mac7.aging <- FindMarkers(Ov.integrated, ident.1 = "Mac7_months3", ident.2 = "Mac7_months15", verbose = FALSE)
head(Mac7.aging, n = 15)
write.csv(Mac7.aging, "Mac7_aging.csv")

T8.aging <- FindMarkers(Ov.integrated, ident.1 = "T8_months3", ident.2 = "T8_months15", verbose = FALSE)
head(T8.aging, n = 15)
write.csv(T8.aging, "T8_aging.csv")

##############################
### AddModuleScore ####
##############################
Ov.integrated = AddModuleScore(object = Ov.integrated,
               features = list(c('Cxcl9', 'Tnfsf10', 'Ccl11', 'Cxcl1', 'Ifng' ,
                                 grep('Tnfa', Ov.integrated@assays$RNA %>% rownames(), value = T)
                                 )),
               name = 'iage_all')

Ov.integrated$iage_all1
print(VlnPlot(object = Ov.integrated, features = 'iage_all1',split.by = "age",group.by = "celltype", 
      pt.size = 0.3))
      
##############################
### AddModuleScore ####
##############################

############################################################
FeaturePlot(Ov.integrated, features = rownames(tail(T8.aging[order(T8.aging$avg_log2FC),],3)), split.by = "age", max.cutoff = 3,
            cols = c("grey", "red"))

sub.cluster.names = Ov.integrated@active.ident[grep('Mac0', Ov.integrated@active.ident)] %>% unique()
sub.cluster.names = factor(sub.cluster.names, levels = sub.cluster.names)
plots <- VlnPlot(Ov.integrated, features = 'Chil3', group.by = 'celltype',split.by = "age",
                 idents = sub.cluster.names,
                 pt.size = 0.3, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(Ov.integrated, features = rownames(tail(T8.aging[order(T8.aging$avg_log2FC),],3)), split.by = "age", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(Ov.integrated, features = c('Cxcl9', 'Tnfsf10', 'Ifng', 'Ccl11', 'Cxcl1'),
# plots <- VlnPlot(Ov.integrated, features = c('Cxcl9', 'Tnfsf10', 'Ifng', 'Ccl11', 'Cxcl1', 
                                             #grep('Tnfa', Ov.integrated@assays$RNA %>% rownames(), value = T)[1:3]),
                 split.by = "age", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(Ov.integrated, features = rownames(tail(Mac7.aging[order(Mac7.aging$avg_log2FC),], 3)), split.by = "age", max.cutoff = 3,
            cols = c("grey", "red"))
plots <- VlnPlot(Ov.integrated, features = rownames(head(T8.aging[order(T8.aging$avg_log2FC),],3)), split.by = "age", group.by = "celltype",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)


###
FeaturePlot(Ov.integrated, features = rownames(tail(Mac0.aging[order(Mac0.aging$avg_log2FC),],3)), split.by = "age", max.cutoff = 3,
            cols = c("grey", "red"))


plots <- VlnPlot(Ov.integrated, features = rownames(tail(Mac0.aging[order(Mac0.aging$avg_log2FC),],3)), split.by = "age", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(Ov.integrated, features = rownames(head(Mac0.aging[order(Mac0.aging$avg_log2FC),],3)), split.by = "age", max.cutoff = 3,
            cols = c("grey", "red"))

plots <- VlnPlot(Ov.integrated, features = rownames(head(Mac0.aging[order(Mac0.aging$avg_log2FC),],3)), split.by = "age", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

###
FeaturePlot(Ov.integrated, features = c('S100a8'), split.by = "age", max.cutoff = 3,
            cols = c("grey", "red"))
############################################################



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
theme_set(theme_cowplot())
Mac2 <- subset(Ov.integrated, idents = "Mac2")
Idents(Mac2) <- "age"
avg.Mac2 <- log1p(AverageExpression(Mac2, verbose = FALSE)$RNA)
avg.Mac2 <- as.data.frame(avg.Mac2)
avg.Mac2$gene <- rownames(avg.Mac2)


cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim"
avg.cd14.mono <- log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA)
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label = c("Cd5l", "Ccl8", "Retnla", "S100a8", "Ifitm1", "Inhba", "Mgl2", "Vcan", "Cxcl3")
p1 <- ggplot(avg.Mac2,aes(avg.Mac2[,1], avg.Mac2[,2])) + geom_point() + ggtitle("Mac2")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
plot_grid(p1, p2)


## ----pressure, echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
plot(pressure)


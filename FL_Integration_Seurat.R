#Author: Faraz Khan, PhD
#FL Integration Seurat

library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)
library(tidyverse)
library(purrr)
library(future)
plan()
plan("multiprocess", workers = 4)
plan()
options(future.globals.maxSize = 7000 * 1024^2)


#Load GEX counts
counts1 <- read.csv('counts.csv', sep = ',', header = T, row.names = 1)
counts1[1:5,1:5]


#Set up Seurat Object
all_merged_cells_obj <- CreateSeuratObject(counts = counts1, min.cells = 3, min.features = 200)
all_merged_cells_obj <- PercentageFeatureSet(object = all_merged_cells_obj, pattern = "^MT-", col.name = "percent.mt")
all_merged_cells_obj
head(all_merged_cells_obj@meta.data)


#Set high bound threshold for nFeature. Data that lies outside 3 times the standard deviation of the median of the data. 
med_T1 <- median(all_merged_cells_obj@meta.data$nFeature_RNA)
stdev_T1 <- sd(all_merged_cells_obj@meta.data$nFeature_RNA)
med_thres_T1 <- med_T1+(stdev_T1*3)


#Subset data
all_merged_cells_obj <- subset(x = all_merged_cells_obj, subset = nFeature_RNA > 200 & nFeature_RNA < med_thres_T1 & nCount_RNA > 200 & percent.mt < 10)
all_merged_cells_obj


#Split the objects to run SCTransform on each sample iteratively
all_merged_cells_obj.list <- SplitObject(all_merged_cells_obj, split.by = "orig.ident")
all_merged_cells_obj.list
for (i in 1:length(all_merged_cells_obj.list)) {
    all_merged_cells_obj.list[[i]] <- SCTransform(all_merged_cells_obj.list[[i]], vars.to.regress = "percent.mt", verbose = T,variable.features.n = 5000)
}
all_merged_cells_obj.list


#Next, select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.
all_merged_cells_obj.features <- SelectIntegrationFeatures(object.list = all_merged_cells_obj.list, nfeatures = 5000)
all_merged_cells_obj.list <- PrepSCTIntegration(object.list = all_merged_cells_obj.list, anchor.features = all_merged_cells_obj.features, verbose = TRUE)
all_merged_cells_obj.list <- lapply(X = all_merged_cells_obj.list, FUN = RunPCA, verbose = FALSE, features = all_merged_cells_obj.features)


#Next, identify anchors and integrate the datasets. Commands are identical to the standard workflow, but make sure to set normalization.method to SCT
all_merged_cells_obj.anchors <- FindIntegrationAnchors(object.list = all_merged_cells_obj.list, normalization.method = "SCT", anchor.features = all_merged_cells_obj.features, reduction = "rpca", verbose = TRUE)
to_integrate <- Reduce(intersect, lapply(all_merged_cells_obj.anchors@object.list, rownames))
all_merged_cells_obj.integrated <- IntegrateData(anchorset = all_merged_cells_obj.anchors, normalization.method = "SCT", verbose = TRUE, features.to.integrate = to_integrate)
all_merged_cells_obj.integrated


#Commands are identical to the standard workflow, but do not run the ScaleData function after integration.

#Remove unwanted genes from variable features
hvgs <- as.data.frame(all_merged_cells_obj.integrated@assays$integrated@var.features)
rownames(hvgs) <- hvgs$`all_merged_cells_obj.integrated@assays$integrated@var.features`
colnames(hvgs) <- 'hvgs'
unwantedg <- as.character(hvgs[rownames(hvgs) %in% c('IGK','IGL','TCRA-VDJ','TCRB-VDJ','TCRG-VDJ','TCRD-VDJ','MT-ATP6','MT-ATP8','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ND1','MT-ND3','MT-ND4','MT-ND4L','MT-ND6','MT-ND2','MT-ND5','MTRNR2L11','MTRNR2L12','MTRNR2L13','MTRNR2L6','MTRNR2L7','MTRNR2L5','MTRNR2L8','MTRNR2L4','MTRNR2L1','MTRNR2L3','MTRNR2L10'),])
wanted_features <- hvgs[!rownames(hvgs) %in% unwantedg,]
wanted_features <- as.data.frame(wanted_features)
all_merged_cells_obj.integrated@assays$integrated@var.features <- as.character(wanted_features$wanted_features)
all_merged_cells_obj.integrated


#Run UMAP, Find Neighbors & Clusters
all_merged_cells_obj.integrated <- RunPCA(all_merged_cells_obj.integrated, npcs = 50, verbose = FALSE)
all_merged_cells_obj.integrated <- RunUMAP(all_merged_cells_obj.integrated, reduction = "pca",dims = 1:18,seed.use = 42)
all_merged_cells_obj.integrated <- FindNeighbors(all_merged_cells_obj.integrated, reduction = "pca",dims = 1:18) 
all_merged_cells_obj.integrated <- FindClusters(all_merged_cells_obj.integrated, resolution = 0.6,random.seed =  42)


#Visualise
DimPlot(object = all_merged_cells_obj.integrated, label = T, pt.size = 0.001, label.size = 6, order = T) + NoLegend()
DimPlot(object = all_merged_cells_obj.integrated, label = F, pt.size = 0.001, group.by = "orig.ident", order = T)


#Switch to RNA assay prior to DE and visualisation
DefaultAssay(all_merged_cells_obj.integrated) <- "RNA"


#Featureplots
fp1<- c("PTPRC", "CD79A", "CD79B", "CD19", "MME")
FeaturePlot(object = all_merged_cells_obj.integrated, features = fp1, pt.size = 1.2, reduction = 'umap',min.cutoff = "q10", max.cutoff = "q90", order = TRUE)


#Run Normalise on RNA assay
all_merged_cells_obj.integrated <- NormalizeData(object = all_merged_cells_obj.integrated, assay="RNA")


#DEA
all_merged_cells_obj.integrated <- AddMetaData(object = all_merged_cells_obj.integrated, as.numeric(all_merged_cells_obj.integrated@meta.data$SeqBatch), col.name = "nSeqBatch") #batch must be numeric
all_markers <- FindAllMarkers(object = all_merged_cells_obj.integrated,only.pos = F,min.pct = 0.25, test.use = "MAST", return.thresh = 0.05, random.seed =  42,assay="RNA",latent.vars="nSeqBatch")
all_markers_adjpval <- all_markers[all_markers$p_val_adj<0.05,]

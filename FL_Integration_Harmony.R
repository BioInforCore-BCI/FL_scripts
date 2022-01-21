#Author: Faraz Khan, PhD
#FL Integration Harmony

library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(data.table)
library(harmony)
library(stringr)
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
head(all_merged_cells_obj@meta.data)
all_merged_cells_obj
all_merged_cells_obj@meta.data


#Set high bound threshold for nFeature. Data that lies outside 3 times the standard deviation of the median of the data. 
med_T1 <- median(all_merged_cells_obj@meta.data$nFeature_RNA)
stdev_T1 <- sd(all_merged_cells_obj@meta.data$nFeature_RNA)
med_thres_T1 <- med_T1+(stdev_T1*3)


#Subset data
all_merged_cells_obj <- subset(x = all_merged_cells_obj, subset = nFeature_RNA > 200 & nFeature_RNA < med_thres_T1 & nCount_RNA > 200 & percent.mt < 10)
all_merged_cells_obj


#Run Normalisation, high variable features and scaling
all_merged_cells_obj <- SCTransform(object = all_merged_cells_obj, vars.to.regress = "percent.mt", verbose = T,variable.features.n = 5000)
all_merged_cells_obj


#Remove unwanted genes from HVG
hvgs <- as.data.frame(all_merged_cells_obj@assays$SCT@var.features)
rownames(hvgs) <- hvgs$`all_merged_cells_obj@assays$SCT@var.features`
colnames(hvgs) <- 'hvgs'
unwantedg <- as.character(hvgs[rownames(hvgs) %in% c('IGK','IGL','TCRA-VDJ','TCRB-VDJ','TCRG-VDJ','TCRD-VDJ','MT-ATP6','MT-ATP8','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ND1','MT-ND3','MT-ND4','MT-ND4L','MT-ND6','MT-ND2','MT-ND5','MTRNR2L11','MTRNR2L12','MTRNR2L13','MTRNR2L6','MTRNR2L7','MTRNR2L5','MTRNR2L8','MTRNR2L4','MTRNR2L1','MTRNR2L3','MTRNR2L10'),])
wanted_features <- hvgs[!rownames(hvgs) %in% unwantedg,]
wanted_features <- as.data.frame(wanted_features)
all_merged_cells_obj@assays$SCT@var.features <- as.character(wanted_features$wanted_features)
all_merged_cells_obj


#Add meta data to Seurat Object
library(stringr)
all_samples <- as.data.frame(rownames(all_merged_cells_obj@meta.data))
colnames(all_samples) <- 'x'
all_samples2 <- str_split(all_samples$x, "_", simplify = T)
all_samples2 <- as.data.frame(all_samples2)
all_merged_cells_obj <- AddMetaData(object = all_merged_cells_obj, all_samples2$V1, col.name = "SampleKey")
all_merged_cells_obj <- AddMetaData(object = all_merged_cells_obj, all_samples2$V2, col.name = "PatientID")
all_merged_cells_obj <- AddMetaData(object = all_merged_cells_obj, all_samples2$V3, col.name = "SampleType")
all_merged_cells_obj <- AddMetaData(object = all_merged_cells_obj, all_samples2$V4, col.name = "CellType")
all_merged_cells_obj <- AddMetaData(object = all_merged_cells_obj, all_samples2$V5, col.name = "SeqBatch")
head(all_merged_cells_obj@meta.data)
tail(all_merged_cells_obj@meta.data)


#Run PCA
all_merged_cells_obj <- RunPCA(all_merged_cells_obj, npcs = 50, verbose = FALSE)


#Batch correction using SeqBatch and PatientID as covariates with Harmony
set.seed(42)
all_merged_cells_obj_har <- RunHarmony(all_merged_cells_obj,c("SeqBatch","PatientID"),assay.use ="SCT", plot_convergence=T, theta=c(2,4),max.iter.cluster=200,max.iter.harmony=100,block.size = 0.01,lambda=c(1,1)) 


#Run UMAP, Find Neighbors & Clusters
all_merged_cells_obj_har <- RunUMAP(all_merged_cells_obj_har, reduction = "harmony",dims = 1:ncol(all_merged_cells_obj_har[["harmony"]]), seed.use = 42)
all_merged_cells_obj_har <- FindNeighbors(all_merged_cells_obj_har, reduction = "harmony",dims = 1:ncol(all_merged_cells_obj_har[["harmony"]])) 
all_merged_cells_obj_har <- FindClusters(all_merged_cells_obj_har, resolution = 1,random.seed =  42)


#Visualise
DimPlot(object = all_merged_cells_obj_har, label = T, pt.size = 0.001, label.size = 3, order = T) + NoLegend()
DimPlot(object = all_merged_cells_obj_har, label = F, pt.size = 0.001, group.by = "SampleType", order = T)
DimPlot(object = all_merged_cells_obj_har, label = F, pt.size = 0.001, group.by = "SampleKey", order = T)
DimPlot(object = all_merged_cells_obj_har, label = F, pt.size = 0.001, group.by = "PatientID", order = T)
DimPlot(object = all_merged_cells_obj_har, label = F, pt.size = 0.001, group.by = "CellType", order = T)


#Featureplots
fp1<- c('PTPRC','CD79A','CD79B','CD3E','CD3D','CD3G','TRDC')
FeaturePlot(object = all_merged_cells_obj_har, features = fp1, pt.size = 1.2, reduction = 'umap',min.cutoff = "q10", max.cutoff = "q90", order = TRUE)


#DEA
all_merged_cells_obj_har <- AddMetaData(object = all_merged_cells_obj_har, as.numeric(all_merged_cells_obj_har@meta.data$SeqBatch), col.name = "nSeqBatch")
all_markers <- FindAllMarkers(object = all_merged_cells_obj_har,only.pos = F,min.pct = 0.25, test.use = "MAST", return.thresh = 0.05,random.seed = 42,latent.vars="nSeqBatch")
all_markers_adjpval <- all_markers[all_markers$p_val_adj<0.05,]

#Author: Faraz Khan, PhD
#FL Integration Harmony subclustering

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


#Load Seurat Object containing all the data
load("SeuratObject.RData")
all_merged_cells_obj_har


#Subset to desired clusters. Idents depends on the user desired cluster numbers
Tcellsobj <- subset(all_merged_cells_obj_har, idents = c("2","3","4","5","9","11","12") , invert=F)


#Re-calculate percent.mt. Note: We use SCT assay for subclustering ONLY. For general pipeline, assay should be "RNA"
Tcellsobj <- PercentageFeatureSet(object = Tcellsobj, pattern = "^MT-", col.name = "percent.mt", assay="SCT")
Tcellsobj
head(Tcellsobj@meta.data)
tail(Tcellsobj@meta.data)


#Subset on percent.mt
Tcellsobj <- subset(x = Tcellsobj, subset = percent.mt < 10)
Tcellsobj


#Run SCTransform
all_merged_cells_obj <- SCTransform(object = Tcellsobj, vars.to.regress = "percent.mt", verbose = T,variable.features.n = 5000, assay="RNA", new.assay.name="SCT")
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


#Run PCA
all_merged_cells_obj <- RunPCA(all_merged_cells_obj, npcs = 20, verbose = FALSE)


#Batch correction using SeqBatch and PatientID as covariates with Harmony
all_merged_cells_obj_har <- RunHarmony(all_merged_cells_obj,c("SeqBatch","PatientID"),assay.use ="SCT", plot_convergence=T, theta=c(2,4),max.iter.cluster=200,max.iter.harmony=100,block.size = 0.01,lambda=c(1,1))


#Run UMAP, Find Neighbors & Clusters
all_merged_cells_obj_har <- RunUMAP(all_merged_cells_obj_har, reduction = "harmony",dims = 1:ncol(all_merged_cells_obj_har[["harmony"]]), seed.use = 42)
all_merged_cells_obj_har <- FindNeighbors(all_merged_cells_obj_har, reduction = "harmony",dims = 1:ncol(all_merged_cells_obj_har[["harmony"]]))
all_merged_cells_obj_har <- FindClusters(all_merged_cells_obj_har, resolution = 0.8,random.seed =  42)


#Visualise
DimPlot(object = all_merged_cells_obj_har, label = T, pt.size = 0.001, label.size = 3, order = T) 
DimPlot(object = all_merged_cells_obj_har, label = F, pt.size = 0.001, group.by = "PatientID", order = T) 
DimPlot(object = all_merged_cells_obj_har, label = F, pt.size = 0.001, group.by = "SampleKey", order = T) 
DimPlot(object = all_merged_cells_obj_har, label = F, pt.size = 0.001, group.by = "SampleType", order = T) 
DimPlot(object = all_merged_cells_obj_har, label = F, pt.size = 0.001, group.by = "CellType", order = T) 


#Featureplots
fp1<- c("PTPRC", "CD3E","CD3D", "CD3G","CD79A", "CD79B", "TRDC","TRGC1","CD8A")
FeaturePlot(object = all_merged_cells_obj_har, features = fp1, pt.size = 1.2, reduction = 'umap',min.cutoff = "q10", max.cutoff = "q90", order = TRUE)


#DEA
all_merged_cells_obj_har <- AddMetaData(object = all_merged_cells_obj_har, as.numeric(all_merged_cells_obj_har@meta.data$SeqBatch), col.name = "nSeqBatch")
all_markers <- FindAllMarkers(object = all_merged_cells_obj_har,only.pos = F,min.pct = 0.25, test.use = "MAST", return.thresh = 0.05,random.seed = 42,latent.vars="nSeqBatch")
all_markers_adjpval <- all_markers[all_markers$p_val_adj<0.05,]

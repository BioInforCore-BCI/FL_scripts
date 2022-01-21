#Author: Faraz Khan, PhD
#FL Seurat

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)


#Load GEX counts
counts1 <- read.csv('counts.csv', sep = ',', header = T, row.names = 1)
counts1[1:5,1:5]


#check how many genes have at least one transcript in each cell
at_least_one <- apply(counts1, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     col="light blue",
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")


#Set up seurat object
SeuratObj <- CreateSeuratObject(counts = counts1, min.cells = 3, min.features = 200)
SeuratObj <- PercentageFeatureSet(object = SeuratObj, pattern = "^MT-", col.name = "percent.mt")
head(SeuratObj@meta.data)
SeuratObj


#We can visualize gene and molecule counts and plot their relationship:
VlnPlot(object = SeuratObj, features = c("nCount_RNA","nFeature_RNA","percent.mt"),log = T,combine = FALSE)


#Set high bound threshold for nFeature. Data that lies outside 3 times the standard deviation of the median of the data.
med_T1 <- median(SeuratObj@meta.data$nFeature_RNA)
stdev_T1 <- sd(SeuratObj@meta.data$nFeature_RNA)
med_thres_T1 <- med_T1+(stdev_T1*3)
plot(SeuratObj@meta.data$nFeature_RNA)
abline(h=med_thres_T1, col="red")


#Subset data
SeuratObj <- subset(x = SeuratObj, subset = nFeature_RNA > 200 & nFeature_RNA < med_thres_T1 & nCount_RNA > 200 & percent.mt < 10)
SeuratObj


#Run Normalisation, high variable features and scaling
SeuratObj <- SCTransform(object = SeuratObj, vars.to.regress = "percent.mt", verbose = T,variable.features.n = 1500)
SeuratObj


#Remove unwanted genes from HVG
hvgs <- as.data.frame(SeuratObj@assays$SCT@var.features)
rownames(hvgs) <- hvgs$`SeuratObj@assays$SCT@var.features`
colnames(hvgs) <- 'hvgs'
unwantedg <- as.character(hvgs[rownames(hvgs) %in% c('IGK','IGL','TCRA-VDJ','TCRB-VDJ','TCRG-VDJ','TCRD-VDJ','MT-ATP6','MT-ATP8','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ND1','MT-ND3','MT-ND4','MT-ND4L','MT-ND6','MT-ND2','MT-ND5','MTRNR2L11','MTRNR2L12','MTRNR2L13','MTRNR2L6','MTRNR2L7','MTRNR2L5','MTRNR2L8','MTRNR2L4','MTRNR2L1','MTRNR2L3','MTRNR2L10'),])
unwantedg
wanted_features <- hvgs[!rownames(hvgs) %in% unwantedg,]
wanted_features <- as.data.frame(wanted_features)
SeuratObj@assays$SCT@var.features <- as.character(wanted_features$wanted_features)
SeuratObj


#Plot histogram before and after normalisation
library(rafalib)
mypar(2)
hist(Matrix::colSums(counts1),
     breaks = 100,
     col = "light blue",
     main = "Total expression per cell before normalisation",
     xlab = "Sum of expression")

hist(Matrix::colSums(SeuratObj),
     breaks = 100,
     col = "light blue",
     main = "Total expression per cell after normalisation",
     xlab = "Sum of expression")


#Now proceed with downstream analysis (i.e. visualization, clustering) 
SeuratObj <- RunPCA(SeuratObj, npcs = 50, verbose = FALSE)


#ElbowPlot for choosing the right dims.
ElbowPlot(object = SeuratObj)


#Examine and visualize PCA results a few different ways
print(x = SeuratObj[["pca"]], dims = 1:5, nfeatures = 5)


#Visualize top genes associated with reduction components
VizDimLoadings(object = SeuratObj, dims = 1:2, reduction = "pca")


#Heatmap
DimHeatmap(object = SeuratObj, dims = 1:20, cells = 500, balanced = TRUE)


#Run UMAP, Find Neighbors & Clusters
SeuratObj <- RunUMAP(SeuratObj, reduction = "pca", dims = 1:18) 
SeuratObj <- FindNeighbors(SeuratObj, reduction = "pca", dims = 1:18) 
SeuratObj <- FindClusters(SeuratObj, resolution = 0.6) 
DimPlot(object = SeuratObj, label = T, pt.size = 0.9)


#Cell Stats
y <- as.data.frame(table(Idents(object = SeuratObj)))
y <- cbind(prop.table(x = table(Idents(object = SeuratObj)))*100,y[,2, drop=F],'sample_name')
names(y) <- c("Cluster","Prop. of cells","Absolute number of cells", "Sample")


#DEA
all_markers <- FindAllMarkers(object = SeuratObj,only.pos = F,min.pct = 0.25, test.use = "MAST", return.thresh = 0.05)
all_markers_adjpval <- all_markers[all_markers$p_val_adj<0.05,]


#Visualise
#FeaturePlot-Tcells
fp<- c("CD3E", "CD3D", "CD3G","CD79A","CD79B", "CD8A", "CD8B", "CD4", "IFNG", "IFI6", "CCL3", "CCL4", "CCL5", "GZMA","GZMH", "PRF1", "FCGR3A", "NKG7", "IL7R", "CCR7","TCF7", "LEF1", "NR3C1", "HAVCR2", "LAG3", "TOX","CXCL13", "PDCD1", "ICOS", "CD200", "CD40LG", "FOXP3","IL2RA","CD44","TRDC")

#FeaturePlot-Bcells
fp <- c("CD19", "MS4A1","IGL","IGK","IGHA","KMT2D","BCL2","HIST1H4C","B2M","IL4I1","HLA-DPA1","ALDOA","PKM","GAPDH","AICDA","LMO2","MIR155HG","FOXO1", "BCL6","MEF2B","CD72","IFIT1","MYC", "REL", "CD83", "NFKBIA","BIRC5","TYMS","TK1","LDHA","STMN1","TOP2A","TNFRSF14","MKI67","CREBBP","CD79A","CD79B","MME")

#Featureplots-nonBnonT
fp<- c("CD79A", "CD79B","CD3E", "CD3D","CD3G", "SDC1", "GZMB", "GZMK", "SELL", "LEF1", "TCF7", "CD68", "C1QA", "APOE","FCGR3A", "NKG7", "KLRF1", "SCT", "LILRA4", "CD34", "DNTT", "IGLL1", "HAVCR2", "LAG3", "LYZ", "S100A9","S100A8", "VCAN", "FCN1", "PTPRC", "CD99", "TRDC", "TRGC1","TNFRSF9", "CD28","COL1A1")

FeaturePlot(object = SeuratObj, features = fp, pt.size = 1.2, reduction = 'umap',min.cutoff = "q10", max.cutoff = "q90", order = TRUE)


#Author: Faraz Khan, PhD
#FL scoring and assignment

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)


#Load the Seurat Object for the data that needs to scored
load("SeuratObject.RData")


#Load the signatures that will be scored. Signatures in columns, genes in rows.
hsig <- read.csv("signatures.csv", check.names = F)
hsig2 <- hsig


#Scale the data
dat <- ScaleData(SeuratObj, features = rownames(SeuratObj), do.scale = F, do.center = T, assay = "RNA")
#Get counts
dat <- as.data.frame(dat@assays$RNA@scale.data)
dat[1:5,1:5]


#Create a new Seurat Object for input in AddModuleScore function
sobj <-  CreateSeuratObject(counts = dat)


#Score
for (i in names(hsig2)){
  p1 = list()
  p1[[i]] <- hsig2[,names(hsig2) == i]
  cd_features <- p1
  sobj <- AddModuleScore(object = sobj,features = cd_features, ctrl = 100, name = i,nbin = 24, seed = 123456, assay = "RNA")
  colnames(sobj@meta.data)[ncol(sobj@meta.data)] <- i
}


#Store the metaprograms scores. Note: Index needs to be adjusted accordingly
head(sobj@meta.data)
mp <- as.data.frame(sobj@meta.data[4:17])


#Make a function for calculating Shannon entropy. Adjust parameters with caution
MP_specificity <- function(x_vector,t=2){
  # Based on Shengli et al 2017
  x <- ifelse(x_vector < 0, 0.01,x_vector+0.01)
  pis <- x/sum(x)
  shannon <- sum(pis * log2(pis))
  ss <- log2(length(x_vector)) - (-shannon)
  H1_2 <- sort(pis, TRUE)[1:2]
  delta_pis <- ifelse(H1_2[1] >= H1_2[2]*t & ss >= 0.8,
                      max(x_vector), "no_Spec")
  return(delta_pis)
}


#Run Shannon entropy and get assignments
mp$get_MP <- colnames(mp)[apply(mp,1,which.max)]
mp$assign <- apply(mp[-ncol(mp)], 1, MP_specificity,t=1.2)
mp$get_MP <- ifelse(mp$assign == 'no_Spec', "No Call", mp$get_MP)
table(mp$get_MP)


#Put the cell assignments as metadata into the original Seurat Object
which(!rownames(mp) == rownames(SeuratObj@meta.data))
SeuratObj <- AddMetaData(SeuratObj, metadata = mp$get_MP, col.name = "CellAssignment")
SeuratObj$CellAssignment <- factor(SeuratObj$CellAssignment, levels = c("No Call","MP1","MP2","MP3","MP4","MP5","MP6","MP7","MP8","MP9","MP10","MP11","MP12","MP13","MP14"))
levels(SeuratObj$CellAssignment)
color_list <- read.csv("colPal_NewMps.csv", header = F)
color_list <- as.character(color_list$V1)
color_list
p1 <- DimPlot(object = SeuratObj, label = F, label.size = 0, pt.size = 0.09, group.by = "CellAssignment", cols = color_list,repel = T) + theme(legend.title = element_text(size = 4),legend.text = element_text(size = 4), legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"),axis.text=element_text(size=6),axis.title=element_text(size=6,face="bold"),axis.line = element_line(colour = 'black', size = 0.3),axis.ticks = element_line(size = 0.07), axis.ticks.length = unit(0.05, "cm"))
p1[[1]]$layers[[1]]$aes_params$alpha = ifelse (SeuratObj@meta.data$CellAssignment == "No Call", .4, 1)
p1


#Individual MP plotting
color_list2 <- color_list
color_list2 <- as.data.frame(color_list2)
color_list2$MP <- c("No Call","MP1","MP2","MP3","MP4","MP5","MP6","MP7","MP8","MP9","MP10","MP11","MP12","MP13","MP14")
for (i in c("No Call","MP1","MP2","MP3","MP4","MP5","MP6","MP7","MP8","MP9","MP10","MP11","MP12","MP13","MP14")){
  mindex <- list()
  mindex[[i]] <- which((color_list2$MP == i))
  usecol <- color_list2
  usecol[1][[1]][-mindex[[i]]] <- 'lightgrey'
  p2 <- DimPlot(object = SeuratObj, label = F, label.size = 0, pt.size = 0.9, group.by = "CellAssignment", cols = as.character(usecol[[1]]),repel = T) + theme_void() + theme(legend.position="none") + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +labs(x = NULL, y = NULL)
  p2[[1]]$layers[[1]]$aes_params$alpha = ifelse (SeuratObj@meta.data$CellAssignment == i, 5, .4)
  pdf(paste(i,"_diagnosticFL.pdf", sep = ""),height = 8, width = 12)
  print(p2)
  dev.off()
}

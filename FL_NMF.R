#Author: Faraz Khan, PhD
#FL NMF

library(Seurat)
library(NMF)


#NMF function
#NMF function written by: Anthony Anene, PhD
F_NMF <- function(data, prefix="NMF", 
                  cluster = 3, nrun=100, 
                  norm=F, ncores=8, 
                  algorithm="brunet", 
                  mode="final", seed=123566){
  
  if (mode == "rank") {
    r <- 2:cluster
    res.r <- NMF::nmf(data, r, algorithm, .opt = paste0("vp", ncores), 
                      nrun = nrun, seed = seed, maxIter = 5000)
    res <- res.r
  } else if (mode == "final") {
    res <- NMF::nmf(data, cluster, algorithm, .opt = paste0("vtp", ncores),
                    nrun = nrun, seed = seed, maxIter = 5000)
  }
  
  return(res)
  # Usage Examples
  # Rank search 
  # rank <- F_NMF(mat, cluster = 5, nrun = 10, norm = F, mode = "rank")
  # Final 
  # Final <- F_NMF(mat, cluster = 3, nrun = 100, norm = F, mode = "final")
}


#Run NMF on each sample normalised counts using only counts from Highly variable gene list (5000)

#Load in counts file

#Run NMF rank survey on counts file to determine optimal number of factors
rank <- F_NMF(ncounts_hvg, cluster = 10, nrun = 10, norm = F, mode = "rank", ncores = 8)


#We can either choose an arbitrary number or find the optimal number from the rank code above
#Run NMF with with 10 factors (arbitrary)
Final <- F_NMF(ncounts_hvg, cluster = 10, nrun = 100, norm = F, mode = "final",ncores = 16)


#Predict and pick top 50 genes from each NMF program
top.n = 50
x <- as.data.frame(NMF::predict(Final, what="features"))
x <- cbind(x, Final@fit@W)
names(x)[1] <- 'features'
#Fold change filter function (Not applied)
diff.vec <- function(y, f=2) {
  y1 <- y[order(y, decreasing=TRUE)][1:2]
  return(ifelse(y1[1] >= y1[2]*f, 1, 0))
}
#Pick top 50 genes per NMF program. Fold change filter is not applied.
x1.top <- lapply(unique(x$features), function(z){
  x1 <- x[x$features == z, ][!names(x) %in% "features"]
  x1 <- x1[order(x1[[which(names(x1) == z)]],
                 decreasing=TRUE), ]
  x1$keep <- apply(x1, 1, FUN = diff.vec,f = 1)
  x1$features <- z
  x1 <- x1[x1$keep == 1, ][1:min(nrow(x1), top.n),]
  return(x1)
})
#get top genes
x1.top <- do.call(rbind,x1.top)
x1.top$gene <- rownames(x1.top)


#Create a combined list of NMF programs generated across all samples. Format: NMF program in column and its genes in rows.


#Score each NMF program using AddModuleScore on the desired Seurat Object
#Load in combined NMF programs list file. Format: NMF program in column and its genes in rows.
alltop <- read.csv("NMF_programs.csv")
#Load in the Seurat Object to score NMF programs on.
load("SeuratObject.RData")
#Run a loop for each unique NMF program
for (i in names(alltop)){
  p1 = list()
  p1[[i]] <- alltop[,names(alltop) == i]
  cd_features <- p1
  all_merged_cells_obj_har <- AddModuleScore(object = all_merged_cells_obj_har,features = cd_features, ctrl = 100, name = i,nbin = 24, seed = 123456)
  colnames(all_merged_cells_obj_har@meta.data)[ncol(all_merged_cells_obj_har@meta.data)] <- i
}
head(all_merged_cells_obj_har@meta.data)


#Correlate NMF programs
library(gplots)
library(amap)
require("RColorBrewer")
myCol <- colorRampPalette(c("blue", "white", "red"))(80)
#Create a data.frame of all the scores across all NMF programs. Index may need to be adjusted according to the user meta data
hc <- as.data.frame(all_merged_cells_obj_har@meta.data[,c(15:114)])
cor1 <- cor(hc)
heatmap.2(cor1,col=myCol,keysize=1.0,key=T,dendrogram="col",margins =c(20,30),lhei = c(0.5,7),lwid = c(0.5,7),scale = "none",notecol="black",density.info="none",trace="none",cexRow=1.3,cexCol=1.6,distfun=function(x) Dist(x, method="pearson"),hclustfun=function(x) hclust(x, method="ward.D2"),reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean))


#Merge NMF programs into NMF Meta Programs based on the correlation cut-off (>=0.5)


#Create a final top30 list of each NMF Meta Program based on the top average NMF weights across correlated NMF programs.


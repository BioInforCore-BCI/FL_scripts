#Author: Faraz Khan, PhD
#FL inferCNV

library(Seurat)
library(infercnv)

##Prepare Annotation file for inferCNV

#Load tumour Bcells sample pre-decontX Seurat Object to extract cell meta data
load('Bcells_SeuratObject.RData')
#Get Idents and cell information. Index may need to be adjusted according to the user meta data
bcells_map <-  as.data.frame(Bcells_SeuratObj@meta.data[,8], row.names = rownames(Bcells_SeuratObj@meta.data))
bcells_map$x <- rownames(bcells_map)
bcells_map <- bcells_map[,c(2,1)]
names(bcells_map) <- c('x','type')
write.table(bcells_map,file="map.txt", sep = "\t", row.names = F, quote = F, col.names = F)


#Load tumour Tcells sample pre-decontX Seurat Object to extract cell meta data
load("Tcells_SeuratObject.RData")
#Get tumour Tcells cell meta data
Tcells_ref_cells <- as.data.frame(colnames(Tcells_SeuratObj),row.names = colnames(Tcells_SeuratObj))
colnames(Tcells_ref_cells) <- 'x'
Tcells_ref_cells['type'] <- 'Tcells'


#Get RLN Bcells cell meta data
RLNBcells <- read.delim("RLN_singletype_cells_map_allsample.txt", header = F) 
names(RLNBcells)[1] <- 'x'
names(RLNBcells)[2] <- 'type'
RLNBcells['type'] <- 'RLN_Bcells'
rownames(RLNBcells) <- RLNBcells$x


#Combine tumour Bcells, tumour Tcells and RLN Bcells meta data to create a Annotation file for input in inferCNV
all_clus_ref <- rbind(Tcells_ref_cells,RLNBcells)
all_clus_ref2 <- rbind(bcells_map,all_clus_ref)
#Save annotation file
write.table(all_clus_ref2,file="annot_clus_RLNadded.txt", sep = "\t", row.names = F, quote = F, col.names = F)



##Load pre-decontx counts file for tumour Bcells, tumour Tcells and RLN Bcells for input in inferCNV

#Load Tumour Bcell patient sample counts file 
counts_TumourBcells <- read.csv('TumourBcells_sample1.csv', sep = ',', header = T, row.names = 1)
#Filter to only tumour Bcells that were used for Seurat clustering and have passed QC.
cell_map <- read.delim("map.txt", header = F, row.names = 1)
counts_TumourBcells <- counts_TumourBcells[, colnames(counts_TumourBcells) %in% rownames(cell_map)]


#Load Tumour Tcell from the same patient counts file that will be used as reference
counts_TumourTcells <- read.csv('TumourTcells_sample1.csv', sep = ',', header = T, row.names = 1)


#Load RLN counts that will be used as reference along with tumour Tcells
counts_RLN <- read.csv('allmergedNLN_dfall_ec.csv', sep = ',', header = T, row.names = 1)
#Load in the file that has RLN Bcell names
RLN_cell_map <- read.delim("RLN_Bcells.txt", header = F, row.names = 1)
#Subset RLN counts file to RLN Bcells.
counts_RLN <- counts_RLN[, colnames(counts_RLN) %in% rownames(RLN_cell_map)]


#Merge all counts
merge.all <- function(x, ..., by = "row.names", all=TRUE) {
  L <- list(...)
  for (i in seq_along(L)) {
    x <- merge(x, L[[i]], by = by, all=T)
    rownames(x) <- x$Row.names
    x$Row.names <- NULL
  }
  return(x)
}
all_merged <- merge.all(counts_TumourBcells,counts_TumourTcells,counts_RLN)
all_merged[is.na(all_merged)] <- 0
all_merged[1:5,1:5]


#Create Seurat Object. Note: This step is optional as inferCNV will perform its own filtering.
seuset <- CreateSeuratObject(
  counts = all_merged,
  min.cells = 3, 
  min.features = 200
)
#Extract counts from Seurat Object and feed in to inferCNV input
all_merged_raw_counts <- as.data.frame(as.matrix(seuset@assays$RNA@counts))
dim(all_merged_raw_counts)


##Run inferCNV (supervised)

#Load in gene order file. This has to match with what the version of annotation used in CellRanger
gene_order <- read.delim("hg38_gene_order_latest_4.txt", sep = '\t', header = F, row.names = 1)
#Load in counts
sc_count_matrix <- all_merged_raw_counts
#Load in Annotations file
cellAnnotations <- read.delim("annot_clus_RLNadded.txt", sep = '\t', header = F, row.names = 1)


#Create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=sc_count_matrix,
                                    annotations_file=cellAnnotations,
                                    delim="\t",
                                    gene_order_file=gene_order,
                                    ref_group_names=c("Tcells","RLN_Bcells"))


# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="sample1_gbc",
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,
                             plot_steps=F,
                             num_threads = 4,
                             HMM_type = "i6",
                             scale_data=F
)



library(Seurat)
library(SeuratData)

setwd("")

obj_dir <- c("sample1","sample2","...")

names(obj_dir) =  c("sample1","sample2","...")


counts <- Read10X(data.dir = obj_dir)

seurat_obj<- CreateSeuratObject(counts,
                                min.feature = 200,
                                min.cells = 3 )


saveRDS(seurat_obj,"scRNA.rds")





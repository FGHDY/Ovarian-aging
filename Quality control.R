library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(batchelor)
library(SeuratWrappers)
library(magrittr)
library(tidyverse)

setwd()

scRNA <- readRDS("scRNA.rds")

#Remove double cells

obj = SplitObject(scRNA, split.by = "orig.ident")
obj_rm=list()
doublets_plot = list()
pc.num = 1:30
dir.create("SingleCell_QC")


RemoveDoublets <-function(
    object,
    doublet.rate,
    pN=0.25,
    pc.num=1:30
){
  sweep.res.list <- paramSweep(object, PCs = pc.num, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  homotypic.prop <- modelHomotypic(object$seurat_clusters)   
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seu.scored <- doubletFinder(object, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                              nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
 
  cname <-colnames(seu.scored[[]])
  DF<-cname[grep('^DF',cname)]
  seu.scored[["doublet"]] <- as.numeric(seu.scored[[DF]]=="Doublet")
  
  # 去除双细胞
  seu.removed <- subset(seu.scored, subset = doublet != 1)
  p1 <- DimPlot(seu.scored, group.by = DF)
  res.list <- list("plot"=p1, "obj"=seu.removed)
  return(res.list)
}
#对每一个样本进行标准聚类的操作
for( i in names(obj)){
  obj[[i]] <- NormalizeData(obj[[i]])
  obj[[i]] <- FindVariableFeatures(obj[[i]], selection.method = "vst", nfeatures = 2000)
  obj[[i]] <- ScaleData(obj[[i]])
  obj[[i]] <- RunPCA(obj[[i]])
  obj[[i]] <- RunUMAP(obj[[i]], dims = 1:30)
  obj[[i]] <- FindNeighbors(obj[[i]], dims = pc.num) %>% FindClusters(resolution = 0.3)
  tmp <- RemoveDoublets(obj[[i]], doublet.rate=0.04,pc.num=pc.num)
  doublets_plot[[i]] <- tmp$plot
}

#Calculate mitochondrial genes
scRNA[["percent.mt"]] = PercentageFeatureSet(scRNA,
                                             pattern = "^MT-")# 人类线粒体基因MT-开头。

#Calculate red blood cell genes
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)]
scRNA[["percent.HB"]] <- PercentageFeatureSet(scRNA,features = HB.genes)

beforeQC_vlnPlot = VlnPlot(scRNA,
                           features = c("nFeature_RNA",
                                        "nCount_RNA",
                                        "percent.mt",
                                        "percent.HB"),
                           ncol = 2,
                           pt.size = 0.05)

dir.create("SingleCell_QC")
ggsave("SingleCell_QC/BeforeQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.pdf", plot = beforeQC_vlnPlot,width = 17,height = 10)
ggsave("SingleCell_QC/BeforeQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.png", plot = beforeQC_vlnPlot,width = 17,height = 10)
beforeQC_vlnPlot    


#Correlation analysis
plot1 <- FeatureScatter(scRNA,feature1 = "nCount_RNA",feature2 = "percent.mt")

plot2 <- FeatureScatter(scRNA,feature1 = "nCount_RNA",feature2 = "nFeature_RNA" )

plot3 <- FeatureScatter(scRNA,feature1 = "nCount_RNA",feature2 = "percent.HB")


#Set quality control conditions
minGene=200 #Number of expressed genes
pctMT=10#Mitochondrial gene ratio
maxHB=5#Red blood cell gene ratio

#qc
scRNA = subset(scRNA,
               subset = nFeature_RNA > minGene & percent.mt < pctMT & percent.HB<maxHB)


afterQC_Vlnplot = VlnPlot(scRNA,
                          features = c("nFeature_RNA",
                                       "nCount_RNA",
                                       "percent.mt",
                                       "percent.HB"),
                          ncol = 2,
                          pt.size = 0.1)



afterQC_Vlnplot
ggsave("SingleCell_QC/afterQC_nFeature_nCount_percent.mt_percent.HB_Vlnplot.pdf", plot = afterQC_Vlnplot,width = 15,height = 8)
ggsave("SingleCell_QC/afterQC_nFeature_nCount_percent.mt_percent.HB_Vlnplot.png", plot = afterQC_Vlnplot,width = 15,height = 8)
afterQC_vlnplot
saveRDS(scRNA,"After_QC.rds")






library(Seurat)
library(DoubletFinder)

RDoublet <- function(tmp){
  sweep.res.list <- paramSweep_v3(tmp, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pKopt <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  pKopt <- pKopt[order(pKopt, decreasing = TRUE) ]
  pKopt <- pKopt[1]
  homotypic.prop <- modelHomotypic(tmp$seurat_clusters)
  nExp_poi <- round(0.05*length(colnames(tmp)))  ## Assuming 5% doublet formation rate
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi, reuse.pANN = FALSE)
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",pKopt,nExp_poi, sep="_"))
  return (tmp)
}

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 3000)
  Seurat <- ScaleData(Seurat, verbose = T, vars.to.regress = c('percent.mt', "percent.rb","S.Score","G2M.Score"))
  
  Seurat <- RunPCA(Seurat, npcs = 100)
  Seurat <- FindNeighbors(Seurat, dims = 1:100)
  Seurat <- FindClusters(Seurat, resolution = 1)
  Seurat <- RunUMAP(Seurat, dims = 1:100)
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}

data_dir <- 'C://Users/rodri/Downloads/GSM6735065'
list.files(data_dir)
D59fetal <- Read10X(data.dir = data_dir)
D59fetalS = CreateSeuratObject(counts = D59fetal)
D59fetalS = CreateSeuratObject(counts = D59fetal, min.cells = 200, min.features = 200)

D59fetalS[["percent.rb"]] <- PercentageFeatureSet(D59fetalS, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA') #change to uppercase for human
D59fetalS[["percent.mt"]] <- PercentageFeatureSet(D59fetalS, pattern = "^MT-") #change to uppercase for human
D59fetalS <- CellCycleScoring(D59fetalS, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE) #remove 'm.' if operating with human data
VlnPlot(D59fetalS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb'), ncol = 4)
D59fetalS <- subset(D59fetalS, subset = nCount_RNA > 300 & nCount_RNA < 50000 & nFeature_RNA > 300 & nFeature_RNA < 6500 & percent.mt < 30 & percent.rb < 40)

D59fetalS <- ProcessSeu(D59fetalS)
D59fetalS <- RDoublet(D59fetalS)
D59fetalS <- subset(D59fetalS, cells = colnames(D59fetalS )[which(D59fetalS [[]][12] == 'Singlet')])
D59fetalS <- subset(D59fetalS , cells = colnames(D59fetalS )[which(D59fetalS [[]][13] == 'Singlet')])
D59fetalS <- ProcessSeu(D59fetalS)
DimPlot(D59fetalS)
D59fetalS$GSM <-  'GSM6735065'
D59fetalS$GSE <- 'GSE218123'
SaveH5Seurat(D59fetalS, 'C://Users/rodri/Downloads/GSE218123_GSM6735065.h5Seurat', overwrite = TRUE)
saveRDS(D59fetalS, 'C://Users/rodri/Downloads/GSE218123_GSM6735065.rds')

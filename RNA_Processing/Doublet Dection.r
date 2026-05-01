library(tidyverse)
library(Seurat)
library(patchwork)
library(sctransform)
library(glmGamPoi)
library(DoubletFinder)

## read data of 164S !
A164S<- Read10X(data.dir = "164S")

#create seurat object
A164S<- CreateSeuratObject(counts = A164S, project = "164S", min.cells = 3, min.features = 200)
View(A164S@meta.data)

#add variable 
A164S$type= "First arthritis" 

#set percent MT
A164S.mt<- PercentageFeatureSet(A164S, pattern = "^MT-", col.name = "percent.mt")

#plot unfiltered data 
unfiltered.164s<- VlnPlot(A164S.mt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfiltered.164s

#filter cells 
Filtered.cells.164s <- subset(A164S.mt, subset = nFeature_RNA >500 
              & nFeature_RNA <5500 & nCount_RNA>500 & nCount_RNA <25000 & percent.mt <15)

#plot filtered cells!  
Filtered.164s<- VlnPlot(Filtered.cells.164s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot both
unfiltered.164s/Filtered.164s

# Normalize cell using standard seurat pipeline
Normalized.164s<- NormalizeData(Filtered.cells.164s)
Normalized.164s <- FindVariableFeatures(Normalized.164s, selection.method = "vst", nfeatures = 2000)
Normalized.164s <- ScaleData(Normalized.164s)
Normalized.164s <- RunPCA(Normalized.164s, verbose = FALSE)

#PC plot
ElbowPlot(Normalized.164s, ndims = 50,reduction = "pca") +
  ggtitle("Number of principle component")  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#Run umap 
Normalized.164s <- RunUMAP(Normalized.164s, dims = 1:40)
DimPlot(Normalized.164s, reduction = "umap", label = TRUE, label.size = 4) +
  theme_void() +
    NoLegend() 

## pK Identification
doublet.A164S<- paramSweep_v3(Normalized.164s, PCs = 1:40, sct = FALSE)
doublet.stats_A164S <- summarizeSweep(doublet.A164S, GT = FALSE)
bcmvn_A164S <- find.pK(doublet.stats_A164S)

#ggplot 
ggplot(bcmvn_A164S, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal()

pK <- bcmvn_A164S %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- Normalized.164s@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.055*nrow(Normalized.164s@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Normalized.164s <- doubletFinder_v3(Normalized.164s, 
                                     PCs = 1:40, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = FALSE)

# visualize doublets
view(Normalized.164s@meta.data)
DimPlot(Normalized.164s, reduction = 'umap', group.by = "DF.classifications_0.25_0.15_409") +
  theme_minimal()

# Number of doublet 
table(Normalized.164s@meta.data$DF.classifications_0.25_0.15_409)
saveRDS(Normalized.164s, file = "Normalized.164s.seurat.RDS")

#filter out doublet ! 
Normalized.164s.singlet<- Normalized.164s[, Normalized.164s@meta.data[, "DF.classifications_0.25_0.15_409"] == "Singlet"]
view(Normalized.164s.singlet)
DimPlot(Normalized.164s.singlet) +
  theme_minimal()

#save singlet object
saveRDS(Normalized.164s.singlet, file = "Single.164S.RDS")

##read data of 164S2 !
A164S2<- Read10X(data.dir = "164S2")

#create seurat object
A164S2<- CreateSeuratObject(counts = A164S2, project = "164S2", min.cells = 3, min.features = 200)

#add variable 
A164S2$type= "Second arthritis" 

#set percent MT
A164S2.mt<- PercentageFeatureSet(A164S2, pattern = "^MT-", col.name = "percent.mt")

#plot unfiltered data 
unfiltered.164S2<- VlnPlot(A164S2.mt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#filter cells 
Filtered.cells.164S2 <- subset(A164S2.mt, subset = nFeature_RNA >500 
              & nFeature_RNA <6000 & nCount_RNA>500 & nCount_RNA <32000 & percent.mt <15)

#plot filtered cells!  
Filtered.164S2<- VlnPlot(Filtered.cells.164S2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot both
unfiltered.164S2/Filtered.164S2
Normalized.164S2<- NormalizeData(Filtered.cells.164S2)
Normalized.164S2 <- FindVariableFeatures(Normalized.164S2, selection.method = "vst", nfeatures = 2000)
Normalized.164S2 <- ScaleData(Normalized.164S2)
Normalized.164S2 <- RunPCA(Normalized.164S2, verbose = FALSE)

#PC plot
ElbowPlot(Normalized.164S2, ndims = 50,reduction = "pca") +
  ggtitle("Number of principle component")  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#Run umap 
Normalized.164S2 <- RunUMAP(Normalized.164S2, dims = 1:40)
DimPlot(Normalized.164S2, reduction = "umap", label = TRUE, label.size = 4) +
  theme_void() +
    NoLegend() 

#Perform doublet finder
## pK Identification
doublet.A164S2<- paramSweep_v3(Normalized.164S2, PCs = 1:40, sct = FALSE)
doublet.stats_A164S2 <- summarizeSweep(doublet.A164S2, GT = FALSE)
bcmvn_A164S2 <- find.pK(doublet.stats_A164S2)

#ggplot 
ggplot(bcmvn_A164S2, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal()

pK <- bcmvn_A164S2 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- Normalized.164S2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.035*nrow(Normalized.164S2@meta.data))  ## Assuming 3.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Normalized.164S2 <- doubletFinder_v3(Normalized.164S2, 
                                     PCs = 1:40, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = FALSE)

# visualize doublets
view(Normalized.164S2@meta.data)
DimPlot(Normalized.164S2, reduction = 'umap', group.by = "DF.classifications_0.25_0.23_145") +
  theme_minimal()

# Number of doublet 
table(Normalized.164S2@meta.data$DF.classifications_0.25_0.23_145)
VlnPlot(Normalized.164S2, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.23_145", pt.size = 0.1) +
  theme_minimal()

#save this object
saveRDS(Normalized.164S2, file = "Normalized.164S2.seurat.RDS")

#filter out doublet ! 
Normalized.164S2.singlet<- Normalized.164S2[, Normalized.164S2@meta.data[, "DF.classifications_0.25_0.23_145"] == "Singlet"]
view(Normalized.164S2.singlet)
DimPlot(Normalized.164S2.singlet) +
  theme_minimal()

#save singlet object
saveRDS(Normalized.164S2.singlet, file = "Single.164S2.RDS")


##read data of 184S2 !
A184S2<- Read10X(data.dir = "184S2")

#create seurat object
A184S2<- CreateSeuratObject(counts = A184S2, project = "184S2", min.cells = 3, min.features = 200)

#add variable 
A184S2$type= "First arthritis" 

#set percent MT
A184S2.mt<- PercentageFeatureSet(A184S2, pattern = "^MT-", col.name = "percent.mt")

#plot unfiltered data 
unfiltered.184S2<- VlnPlot(A184S2.mt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#filter cells 
Filtered.cells.184S2 <- subset(A184S2.mt, subset = nFeature_RNA >500 
              & nFeature_RNA <6000 & nCount_RNA>500 & nCount_RNA <30000 & percent.mt <15)

#plot filtered cells!  
Filtered.184S2<- VlnPlot(Filtered.cells.184S2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Filtered.184S2

#plot both
unfiltered.184S2/Filtered.184S2

Normalized.184S2<- NormalizeData(Filtered.cells.184S2)
Normalized.184S2 <- FindVariableFeatures(Normalized.184S2, selection.method = "vst", nfeatures = 2000)
Normalized.184S2 <- ScaleData(Normalized.184S2)
Normalized.184S2 <- RunPCA(Normalized.184S2, verbose = FALSE)

#PC plot
ElbowPlot(Normalized.184S2, ndims = 50,reduction = "pca") +
  ggtitle("Number of principle component")  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#Run umap 
Normalized.184S2 <- RunUMAP(Normalized.184S2, dims = 1:40)
DimPlot(Normalized.184S2, reduction = "umap", label = TRUE, label.size = 4) +
  theme_void() +
    NoLegend() 

#Perform doublet finder
## pK Identification
library(DoubletFinder)
doublet.A184S2<- paramSweep_v3(Normalized.184S2, PCs = 1:40, sct = FALSE)
doublet.stats_A184S2 <- summarizeSweep(doublet.A184S2, GT = FALSE)
bcmvn_A184S2 <- find.pK(doublet.stats_A184S2)

#ggplot 
ggplot(bcmvn_A184S2, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal()

pK <- bcmvn_A184S2 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- Normalized.184S2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations 
nExp_poi <- round(0.023*nrow(Normalized.184S2@meta.data))  ## Assuming 2.3 percent doublet
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Normalized.184S2 <- doubletFinder_v3(Normalized.184S2, 
                                     PCs = 1:40, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = FALSE)

# visualize doublets
view(Normalized.184S2@meta.data)
DimPlot(Normalized.184S2, reduction = 'umap', group.by = "DF.classifications_0.25_0.14_82") +
  theme_minimal()

# Number of doublet 
table(Normalized.184S2@meta.data$DF.classifications_0.25_0.14_82)

VlnPlot(Normalized.184S2, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.14_82", pt.size = 0.1) +
  theme_minimal()

#save this object
saveRDS(Normalized.184S2, file = "Normalized.184S2.seurat.RDS")

#filter out doublet ! 
Normalized.184S2.singlet<- Normalized.184S2[, Normalized.184S2@meta.data[, "DF.classifications_0.25_0.14_82"] == "Singlet"]
view(Normalized.184S2.singlet)
DimPlot(Normalized.184S2.singlet) +
  theme_minimal()

#save singlet object
saveRDS(Normalized.184S2.singlet, file = "Single.184S2.RDS")

## read data of 184S3 !
A184S3<- Read10X(data.dir = "184S3")

#create seurat object
A184S3<- CreateSeuratObject(counts = A184S3, project = "184S3", min.cells = 3, min.features = 200)
A184S3

#add variable 
A184S3$type= "Second arthritis" 

#set percent MT
A184S3.mt<- PercentageFeatureSet(A184S3, pattern = "^MT-", col.name = "percent.mt")

#plot unfiltered data 
unfiltered.184S3<- VlnPlot(A184S3.mt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfiltered.184S3

#filter cells 
Filtered.cells.184S3 <- subset(A184S3.mt, subset = nFeature_RNA >500 
              & nFeature_RNA <6000 & nCount_RNA>500 & nCount_RNA <32000 & percent.mt <15)

#plot filtered cells!  
Filtered.184S3<- VlnPlot(Filtered.cells.184S3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Filtered.184S3

#plot both
unfiltered.184S3/Filtered.184S3

Normalized.184S3<- NormalizeData(Filtered.cells.184S3)
Normalized.184S3 <- FindVariableFeatures(Normalized.184S3, selection.method = "vst", nfeatures = 2000)
Normalized.184S3 <- ScaleData(Normalized.184S3)
Normalized.184S3 <- RunPCA(Normalized.184S3, verbose = FALSE)

#PC plot
ElbowPlot(Normalized.184S3, ndims = 50,reduction = "pca") +
  ggtitle("Number of principle component")  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#Run umap 
Normalized.184S3 <- RunUMAP(Normalized.184S3, dims = 1:40)
DimPlot(Normalized.184S3, reduction = "umap", label = TRUE, label.size = 4) +
  theme_void() +
    NoLegend() 

#Perform doublet finder
## pK Identification
doublet.A184S3<- paramSweep_v3(Normalized.184S3, PCs = 1:40, sct = FALSE)
doublet.stats_A184S3 <- summarizeSweep(doublet.A184S3, GT = FALSE)
bcmvn_A184S3 <- find.pK(doublet.stats_A184S3)

#ggplot 
ggplot(bcmvn_A184S3, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal()

pK <- bcmvn_A184S3 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- Normalized.184S3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations 
nExp_poi <- round(0.015*nrow(Normalized.184S3@meta.data))  ## Assuming 2.3 percent doublet
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Normalized.184S3 <- doubletFinder_v3(Normalized.184S3, 
                                     PCs = 1:40, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = FALSE)

# visualize doublets
view(Normalized.184S3@meta.data)
DimPlot(Normalized.184S3, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_27") +
  theme_minimal()

# Number of doublet 
table(Normalized.184S3@meta.data$DF.classifications_0.25_0.3_27)

VlnPlot(Normalized.184S3, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.3_27", pt.size = 0.1) +
  theme_minimal()

#save this object

saveRDS(Normalized.184S3, file = "Normalized.184S3.seurat.RDS")

#filter out doublet ! 
Normalized.184S3.singlet<- Normalized.184S3[, Normalized.184S3@meta.data[, "DF.classifications_0.25_0.3_27"] == "Singlet"]
view(Normalized.184S3.singlet)
DimPlot(Normalized.184S3.singlet) +
  theme_minimal()

#save singlet object
saveRDS(Normalized.184S3.singlet, file = "Single.184S3.RDS")

       
##read data of 218S !
A218S<- Read10X(data.dir = "218S")

#create seurat object
A218S<- CreateSeuratObject(counts = A218S, project = "218S", min.cells = 3, min.features = 200)

#add variable 
A218S$type= "First arthritis" 

#set percent MT
A218S.mt<- PercentageFeatureSet(A218S, pattern = "^MT-", col.name = "percent.mt")

#plot unfiltered data 
unfiltered.218S<- VlnPlot(A218S.mt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfiltered.218S

#filter cells 
Filtered.cells.218S <- subset(A218S.mt, subset = nFeature_RNA >500 
              & nFeature_RNA <5900 & nCount_RNA>500 & nCount_RNA <25000 & percent.mt <15)

#plot filtered cells!  
Filtered.218S<- VlnPlot(Filtered.cells.218S, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Filtered.218S

#plot both
unfiltered.218S/Filtered.218S
Normalized.218S<- NormalizeData(Filtered.cells.218S)
Normalized.218S<- FindVariableFeatures(Normalized.218S, selection.method = "vst", nfeatures = 2000)
Normalized.218S <- ScaleData(Normalized.218S)
Normalized.218S <- RunPCA(Normalized.218S, verbose = FALSE)

#PC plot
ElbowPlot(Normalized.218S, ndims = 50,reduction = "pca") +
  ggtitle("Number of principle component")  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#Run umap 
Normalized.218S <- RunUMAP(Normalized.218S, dims = 1:40)
DimPlot(Normalized.218S, reduction = "umap", label = TRUE, label.size = 4) +
  theme_void() +
    NoLegend() 

#Perform doublet finder
## pK Identification
doublet.A218S<- paramSweep_v3(Normalized.218S, PCs = 1:40, sct = FALSE)
doublet.stats_A218S <- summarizeSweep(doublet.A218S, GT = FALSE)
bcmvn_A218S <- find.pK(doublet.stats_A218S)

#ggplot 
ggplot(bcmvn_A218S, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal()

pK <- bcmvn_A218S %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- Normalized.218S@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations 
nExp_poi <- round(0.06*nrow(Normalized.218S@meta.data))  ## Assuming 2.3 percent doublet
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Normalized.218S <- doubletFinder_v3(Normalized.218S, 
                                     PCs = 1:40, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = FALSE)

# visualize doublets
view(Normalized.218S@meta.data)
DimPlot(Normalized.218S, reduction = 'umap', group.by = "DF.classifications_0.25_0.26_515") +
  theme_minimal()

# Number of doublet 
table(Normalized.218S@meta.data$DF.classifications_0.25_0.26_515)

VlnPlot(Normalized.218S, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.26_515", pt.size = 0.1) +
  theme_minimal()

#save this object
saveRDS(Normalized.218S, file = "Normalized.218S.seurat.RDS")

#filter out doublet ! 
Normalized.218S.singlet<- Normalized.218S[, Normalized.218S@meta.data[, "DF.classifications_0.25_0.26_515"] == "Singlet"]
view(Normalized.218S.singlet)
DimPlot(Normalized.218S.singlet) +
  theme_minimal()

#save singlet object
saveRDS(Normalized.218S.singlet, file = "Single.218S.RDS")
                                       
#read data of 218S2 !
A218S2<- Read10X(data.dir = "218S2")

#create seurat object
A218S2<- CreateSeuratObject(counts = A218S2, project = "218S2", min.cells = 3, min.features = 200)
A218S2

#add variable 
A218S2$type= "Second arthritis" 

#set percent MT
A218S2.mt<- PercentageFeatureSet(A218S2, pattern = "^MT-", col.name = "percent.mt")

#plot unfiltered data 
unfiltered.218S2<- VlnPlot(A218S2.mt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#filter cells 
Filtered.cells.218S2 <- subset(A218S2.mt, subset = nFeature_RNA >500 
              & nFeature_RNA <6000 & nCount_RNA>500 & nCount_RNA <28000 & percent.mt <15)

#plot filtered cells!  
Filtered.218S2<- VlnPlot(Filtered.cells.218S2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Filtered.218S2

#plot both
unfiltered.218S2/Filtered.218S2
Normalized.218S2<- NormalizeData(Filtered.cells.218S2)
Normalized.218S2<- FindVariableFeatures(Normalized.218S2, selection.method = "vst", nfeatures = 2000)
Normalized.218S2 <- ScaleData(Normalized.218S2)
Normalized.218S2 <- RunPCA(Normalized.218S2, verbose = FALSE)

#PC plot
ElbowPlot(Normalized.218S2, ndims = 50,reduction = "pca") +
  ggtitle("Number of principle component")  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#Run umap 
Normalized.218S2 <- RunUMAP(Normalized.218S2, dims = 1:40)
DimPlot(Normalized.218S2, reduction = "umap", label = TRUE, label.size = 4) +
  theme_void() +
    NoLegend() 

#Perform doublet finder
## pK Identification
doublet.A218S2<- paramSweep_v3(Normalized.218S2, PCs = 1:40, sct = FALSE)
doublet.stats_A218S2 <- summarizeSweep(doublet.A218S2, GT = FALSE)
bcmvn_A218S2 <- find.pK(doublet.stats_A218S2)

#ggplot 
ggplot(bcmvn_A218S2, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal()

pK <- bcmvn_A218S2 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- Normalized.218S2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations 
nExp_poi <- round(0.031*nrow(Normalized.218S2@meta.data))  ## Assuming 3.1 percent doublet
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Normalized.218S2 <- doubletFinder_v3(Normalized.218S2, 
                                     PCs = 1:40, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = FALSE)

# visualize doublets
view(Normalized.218S2@meta.data)
DimPlot(Normalized.218S2, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_119") +
  theme_minimal()

# Number of doublet 
table(Normalized.218S2@meta.data$DF.classifications_0.25_0.005_119)

VlnPlot(Normalized.218S2, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.005_119", pt.size = 0.1) +
  theme_minimal()

#save this object
saveRDS(Normalized.218S2, file = "Normalized.218S2.seurat.RDS")

#filter out doublet ! 
Normalized.218S2.singlet<- Normalized.218S2[, Normalized.218S2@meta.data[, "DF.classifications_0.25_0.005_119"] == "Singlet"]
view(Normalized.218S2.singlet)
DimPlot(Normalized.218S2.singlet) +
  theme_minimal()

#save singlet object
saveRDS(Normalized.218S2.singlet, file = "Single.218S2.RDS")

#read data of SA5 !
ASA5<- Read10X(data.dir = "SA5")

#create seurat object
ASA5<- CreateSeuratObject(counts = ASA5, project = "SA5", min.cells = 3, min.features = 200)

#add variable 
ASA5$type= "Osteoarthritis" 

#set percent MT
ASA5.mt<- PercentageFeatureSet(ASA5, pattern = "^MT-", col.name = "percent.mt")

#plot unfiltered data 
unfiltered.SA5<- VlnPlot(ASA5.mt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#filter cells 
Filtered.cells.SA5 <- subset(ASA5.mt, subset = nFeature_RNA >500 
              & nFeature_RNA <6900 & nCount_RNA>500 & nCount_RNA <31000 & percent.mt <15)

#plot filtered cells!  
Filtered.SA5<- VlnPlot(Filtered.cells.SA5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Filtered.SA5

#plot both
unfiltered.SA5/Filtered.SA5
Normalized.SA5<- NormalizeData(Filtered.cells.SA5)
Normalized.SA5<- FindVariableFeatures(Normalized.SA5, selection.method = "vst", nfeatures = 2000)
Normalized.SA5 <- ScaleData(Normalized.SA5)
Normalized.SA5 <- RunPCA(Normalized.SA5, verbose = FALSE)

#PC plot
ElbowPlot(Normalized.SA5, ndims = 50,reduction = "pca") +
  ggtitle("Number of principle component")  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#Run umap 
Normalized.SA5 <- RunUMAP(Normalized.SA5, dims = 1:40)
DimPlot(Normalized.SA5, reduction = "umap", label = TRUE, label.size = 4) +
  theme_void() +
    NoLegend() 

#Perform doublet finder
## pK Identification
doublet.ASA5<- paramSweep_v3(Normalized.SA5, PCs = 1:40, sct = FALSE)
doublet.stats_ASA5 <- summarizeSweep(doublet.ASA5, GT = FALSE)
bcmvn_ASA5 <- find.pK(doublet.stats_ASA5)

#ggplot 
ggplot(bcmvn_ASA5, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal()

pK <- bcmvn_ASA5 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- Normalized.SA5@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations 
nExp_poi <- round(0.054*nrow(Normalized.SA5@meta.data))  ## Assuming 3.1 percent doublet
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Normalized.SA5 <- doubletFinder_v3(Normalized.SA5, 
                                     PCs = 1:40, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = FALSE)

# visualize doublets
view(Normalized.SA5@meta.data)
DimPlot(Normalized.SA5, reduction = 'umap', group.by = "DF.classifications_0.25_0.27_356") +
  theme_minimal()

# Number of doublet 
table(Normalized.SA5@meta.data$DF.classifications_0.25_0.27_356)

VlnPlot(Normalized.SA5, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.27_356", pt.size = 0.1) +
  theme_minimal()

#save this object
saveRDS(Normalized.SA5, file = "Normalized.SA5.seurat.RDS")

#filter out doublet ! 
Normalized.SA5.singlet<- Normalized.SA5[, Normalized.SA5@meta.data[, "DF.classifications_0.25_0.27_356"] == "Singlet"]
view(Normalized.SA5.singlet)
DimPlot(Normalized.SA5.singlet) +
  theme_minimal()


#save singlet object
saveRDS(Normalized.SA5.singlet, file = "Single.SA5.RDS")
                                              
##read data of SA6 !
ASA6<- Read10X(data.dir = "SA6")

#create seurat object
ASA6<- CreateSeuratObject(counts = ASA6, project = "SA6", min.cells = 3, min.features = 200)
ASA6

#add variable 
ASA6$type= "Osteoarthritis" 

#set percent MT
ASA6.mt<- PercentageFeatureSet(ASA6, pattern = "^MT-", col.name = "percent.mt")

#plot unfiltered data 
unfiltered.SA6<- VlnPlot(ASA6.mt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfiltered.SA6

#filter cells 
Filtered.cells.SA6 <- subset(ASA6.mt, subset = nFeature_RNA >500 
              & nFeature_RNA <7500
              & nCount_RNA>500 & nCount_RNA <58000 & percent.mt <15)

#plot filtered cells!  
Filtered.SA6<- VlnPlot(Filtered.cells.SA6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Filtered.SA6

#plot both
unfiltered.SA6/Filtered.SA6

Normalized.SA6<- NormalizeData(Filtered.cells.SA6)
Normalized.SA6<- FindVariableFeatures(Normalized.SA6, selection.method = "vst", nfeatures = 2000)
Normalized.SA6 <- ScaleData(Normalized.SA6)
Normalized.SA6 <- RunPCA(Normalized.SA6, verbose = FALSE)

#PC plot
ElbowPlot(Normalized.SA6, ndims = 50,reduction = "pca") +
  ggtitle("Number of principle component")  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#Run umap 
Normalized.SA6 <- RunUMAP(Normalized.SA6, dims = 1:40)
DimPlot(Normalized.SA6, reduction = "umap", label = TRUE, label.size = 4) +
  theme_void() +
    NoLegend() 

#Perform doublet finder
## pK Identification
doublet.ASA6<- paramSweep_v3(Normalized.SA6, PCs = 1:40, sct = FALSE)
doublet.stats_ASA6 <- summarizeSweep(doublet.ASA6, GT = FALSE)
bcmvn_ASA6 <- find.pK(doublet.stats_ASA6)

#ggplot 
ggplot(bcmvn_ASA6, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal()

pK <- bcmvn_ASA6 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- Normalized.SA6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations 
nExp_poi <- round(0.05*nrow(Normalized.SA6@meta.data))  ## Assuming 3.1 percent doublet
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Normalized.SA6 <- doubletFinder_v3(Normalized.SA6, 
                                     PCs = 1:45, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = FALSE)

# visualize doublets
view(Normalized.SA6@meta.data)
DimPlot(Normalized.SA6, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_307") +
  theme_minimal()

# Number of doublet 
table(Normalized.SA6@meta.data$DF.classifications_0.25_0.3_307)

VlnPlot(Normalized.SA6, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.3_307", pt.size = 0.1) +
  theme_minimal()

#save this object
saveRDS(Normalized.SA6, file = "Normalized.SA6.seurat.RDS")



#filter out doublet ! 
Normalized.SA6.singlet<- Normalized.SA6[, Normalized.SA6@meta.data[, "DF.classifications_0.25_0.3_307"] == "Singlet"]
view(Normalized.SA6.singlet)
DimPlot(Normalized.SA6.singlet) +
  theme_minimal()

#save singlet object
saveRDS(Normalized.SA6.singlet, file = "Single.SA6.RDS")
                                  
#read data of SA6 !
ASA7<- Read10X(data.dir = "SA7")

#create seurat object
ASA7<- CreateSeuratObject(counts = ASA7, project = "SA7", min.cells = 3, min.features = 200)

#add variable 
ASA7$type= "Osteoarthritis" 

#set percent MT
ASA7.mt<- PercentageFeatureSet(ASA7, pattern = "^MT-", col.name = "percent.mt")

#plot unfiltered data 
unfiltered.SA7<- VlnPlot(ASA7.mt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
unfiltered.SA7

#filter cells 
Filtered.cells.SA7 <- subset(ASA7.mt, subset = nFeature_RNA >500 
              & nFeature_RNA <7500
              & nCount_RNA>500 & nCount_RNA <50000 & percent.mt <15)

#plot filtered cells!  
Filtered.SA7<- VlnPlot(Filtered.cells.SA7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot both
unfiltered.SA7/Filtered.SA7

Normalized.SA7<- NormalizeData(Filtered.cells.SA7)
Normalized.SA7<- FindVariableFeatures(Normalized.SA7, selection.method = "vst", nfeatures = 2000)
Normalized.SA7 <- ScaleData(Normalized.SA7)
Normalized.SA7 <- RunPCA(Normalized.SA7, verbose = FALSE)

#PC plot
ElbowPlot(Normalized.SA7, ndims = 50,reduction = "pca") +
  ggtitle("Number of principle component")  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#Run umap 
Normalized.SA7 <- RunUMAP(Normalized.SA7, dims = 1:40)
DimPlot(Normalized.SA7, reduction = "umap", label = TRUE, label.size = 4) +
  theme_void() +
    NoLegend() 

#Perform doublet finder
## pK Identification
doublet.ASA7<- paramSweep_v3(Normalized.SA7, PCs = 1:40, sct = FALSE)
doublet.stats_ASA7 <- summarizeSweep(doublet.ASA7, GT = FALSE)
bcmvn_ASA7 <- find.pK(doublet.stats_ASA7)

#ggplot 
ggplot(bcmvn_ASA7, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal()

pK <- bcmvn_ASA7 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- Normalized.SA7@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations 
nExp_poi <- round(0.05*nrow(Normalized.SA7@meta.data))  ## Assuming 3.1 percent doublet
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Normalized.SA7 <- doubletFinder_v3(Normalized.SA7, 
                                     PCs = 1:40, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = FALSE)

# visualize doublets
view(Normalized.SA7@meta.data)
DimPlot(Normalized.SA7, reduction = 'umap', group.by = "DF.classifications_0.25_0.26_139") +
  theme_minimal()

# Number of doublet 
table(Normalized.SA7@meta.data$DF.classifications_0.25_0.26_139)


VlnPlot(Normalized.SA7, features = "nFeature_RNA", group.by = "DF.classifications_0.25_0.26_139", pt.size = 0.1) +
  theme_minimal()

#save this object
saveRDS(Normalized.SA7, file = "Normalized.SA7.seurat.RDS")


#filter out doublet ! 

Normalized.SA7.singlet<- Normalized.SA7[, Normalized.SA7@meta.data[, "DF.classifications_0.25_0.26_139"] == "Singlet"]
view(Normalized.SA7.singlet)
DimPlot(Normalized.SA7.singlet) +
  theme_minimal()

#save singlet object
saveRDS(Normalized.SA7.singlet, file = "Single.SA7.RDS")








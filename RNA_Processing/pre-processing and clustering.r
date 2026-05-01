library(tidyverse)
library(Seurat)
library(patchwork)
library(sctransform)
library(glmGamPoi)
library(harmony)

# Read Singlet objects 
p_164S<- readRDS("Singlet.164S.RDS")
p_164S2<- readRDS("Singlet.164S2.RDS")
p_184S2<- readRDS("Singlet.184S2.RDS")
p_184S3<- readRDS("Singlet.184S3.RDS")
p_218S<- readRDS("Singlet.218S.RDS")
p_218S2<- readRDS("Singlet.218S2.RDS")
p_SA5<- readRDS("Singlet.SA5.RDS")
p_SA6<- readRDS("Singlet.SA6.RDS")
p_SA7<- readRDS("Singlet.SA7.RDS")

# Merge the singlet objects
Singlet<- merge(Normalized.164s.singlet, y=c(p_164S, p_164S2, p_184S2, p_184S3, p_218S, p_218S2, p_SA5, p_SA6, p_SA7))
# Save merged singlet object! 
saveRDS(Singlet, file = "Singlet_SCT_merge.RDS")

#plot unfiltered data 
VlnPlot(Singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Feature and feature relationship 
plot1 <- FeatureScatter(Singlet, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) 
plot2 <- FeatureScatter(Singlet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
plot1|plot2 

# Normalize the data using SCTransform
#SCTransform
Arthritis<- SCTransform(Singlet, 
                        method= "glmGamPoi", 
                        vars.to.regress = "percent.mt", 
                        verbose= FALSE)

#Run PCA
Arthritis<- RunPCA(Arthritis, verbose = FALSE)

#Save pca object 
saveRDS(Arthritis, file = "preprocess_SCT_PCA.RDS")

#Elbow plot
ElbowPlot(Arthritis, ndims = 50,reduction = "pca") +
  ggtitle("Number of principle component")  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Heatmap for Deg
DimHeatmap(Arthritis, dims = 1:3, cells = 500, balanced = TRUE)

# Integration and clustering 
library(harmony)
Arthritis<- readRDS("preprocess_SCT_PCA.RDS")
Arthritis.integration<- RunHarmony(Arthritis ,group.by.vars = "orig.ident")

#save harmony object
saveRDS(Arthritis.integration, file = "integration.harmony.RDS")

UMAP.singlparam<- FindNeighbors(Arthritis.integration, reduction = "harmony", dims = 1:43, verbose = FALSE)
UMAP.singlparam<- FindClusters(UMAP.singlparam, resolution = 0.4)
UMAP.singlparam<- RunUMAP(object = UMAP.singlparam,
                          reduction = "harmony",
                          dims = 1:43)

#save this object
saveRDS(UMAP.singlparam, file = "UMAP_SCT_Filter.RDS")

#plot umap 
UMAP.singlparam<- readRDS("UMAP_SCT_Filter.RDS")
DimPlot(UMAP.singlparam, reduction= "umap", label = TRUE) +
  theme_minimal() +
  NoLegend()

# Visualize All the three groups in UMAP
DimPlot(UMAP.singlparam, reduction = "umap", split.by = "type") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  NoLegend()

#' T cell subclustering ============================================================================================

#Select T cell cluster*
T.cell.cluster <- subset(seurat.obj, idents = c("0", "2", "3", "4", "7", "9", "11", "12", "13", "14"))

# plot T cell cluster
DimPlot(T.cell.cluster, reduction = "umap", label = FALSE, label.size = 3) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14), 
        axis.title= element_text(size = 14, face = "bold"), 
        axis.text = element_text(size=13, face = "bold"), 
        axis.text.x = element_text(angle = 90, vjust = 0.1)) +
  labs(y= "UMAP__2", x= "UMAP__1") + NoAxes()


# clean up the debris 
T.subcluster.plot<- DimPlot(T.cell.cluster, reduction= "umap", label = TRUE) +
 theme_minimal() +
  NoLegend()
seurat.T.cluster<- CellSelector(plot = T.subcluster.plot, T.cell.cluster)
DimPlot(seurat.T.cluster, reduction= "umap", label = TRUE) +
  theme_minimal() +
  NoLegend()

# Select cell again
T.cell.cluster <- subset(seurat.T.cluster, idents = c("0", "2", "3", "4", "7", "9", "11", "12", "13", "14"))
DimPlot(T.cell.cluster)

T.subcluster.plot<- DimPlot(T.cell.cluster, reduction= "umap", label = TRUE) +
 theme_minimal() +
  NoLegend()
seurat.T.cluster<- CellSelector(plot = T.subcluster.plot, T.cell.cluster)

DimPlot(seurat.T.cluster, reduction= "umap", label = TRUE) +
  theme_minimal() +
  NoLegend()

T.cell.cluster <- subset(seurat.T.cluster, idents = c("0", "2", "3", "4", "7", "9", "11", "12", "13", "14"))
t.sub<- DimPlot(T.cell.cluster, label = TRUE) +
  theme_void() +
  NoLegend() +
  ggtitle("subset all CD3 clusters") +
  theme(plot.title = element_text(hjust = 0.5))
t.sub

# SCTransform
T.subcluster<- SCTransform(T.cell.cluster, 
                         method= "glmGamPoi")

T.subcluster<- RunPCA(T.subcluster, verbose = FALSE)
ElbowPlot(T.subcluster, ndims = 50,reduction = "pca") +
  ggtitle("Number of principle component")  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

T.cell.integration<- RunHarmony(T.subcluster,
                              group.by.vars = "orig.ident")

sub.cluster.T.pc15<- FindNeighbors(T.cell.integration, 
                                 reduction = "harmony", 
                               dims = 1:15)

sub.cluster.T.pc15<- FindClusters(sub.cluster.T.pc15,
                                 resolution = 0.6)

sub.cluster.T.pc15<- RunUMAP(object = sub.cluster.T.pc15,
                          reduction = "harmony",
                          dims = 1:15)

#plot umap 
sub.cluster.T<- read_rds("sub.cluster.T.rds")

#umap plot 
DimPlot(sub.cluster.T.pc15, reduction= "umap", label = TRUE) +
  theme_void() +
  NoLegend() +
  NoAxes() + 
  ggtitle("PC= 15") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

saveRDS(sub.cluster.T, file = "umap.obj.pc15.rds")


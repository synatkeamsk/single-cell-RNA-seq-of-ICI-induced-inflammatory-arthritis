library(tidyverse)
library(CellChat)
library(Seurat)
library(reticulate)
library(CellChat)
library(patchwork)
library(circlize)

seurat.obj.cellchat<- readRDS("all.clusters.annotated.obj.rds")
DimPlot(seurat.obj.cellchat, reduction = "umap", label = TRUE, label.size = 3.5) +
  theme_minimal() +
  NoLegend()

#' Filter first arthritis obj
seurat.obj.first.art<- subset(seurat.obj.cellchat, subset = type == "First arthritis")
View(seurat.obj.first.art@meta.data)

# Create the NL CellChat object 
cellchat.first.art <- createCellChat(object = seurat.obj.first.art, group.by = "ident", assay = "SCT")

# Add the Secreted Signaling database in the CellChat object
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchat.first.art@DB <- CellChatDB.use

# ubset and pre-processing the expression data 
cellchat.first.art <- subsetData(cellchat.first.art)
cellchat.first.art <- identifyOverExpressedGenes(cellchat.first.art)
cellchat.first.art <- identifyOverExpressedInteractions(cellchat.first.art)

# project gene expression data onto protein-protein interaction (PPI)
cellchat.first.art <- projectData(cellchat.first.art, PPI.human) 

# Compute the communication probability and infer cellular communication network
cellchat.first.art <- computeCommunProb(cellchat.first.art, raw.use = FALSE) 

# Filter out the cell-cell communication min.cells = 10
cellchat.first.art <- filterCommunication(cellchat.first.art, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat.first.art <- computeCommunProbPathway(cellchat.first.art)

# Calculate the aggregated cell-cell communication network
cellchat.first.art <- aggregateNet(cellchat.first.art)

# Compute the network centrality scores
cellchat.first.art <- netAnalysis_computeCentrality(cellchat.first.art, slot.name = "netP")

#' ===============================================================================================

#' Filter second arthritis obj
seurat.obj.second.art<- subset(seurat.obj.cellchat, subset = type == "Second arthritis")

#  Create the NL CellChat object 
cellchat.second.art <- createCellChat(object = seurat.obj.second.art, group.by = "ident", assay = "SCT")

# Add the Secreted Signaling database in the CellChat object
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchat.second.art@DB <- CellChatDB.use

# Subset and pre-processing the expression data 
cellchat.second.art <- subsetData(cellchat.second.art)
cellchat.second.art <- identifyOverExpressedGenes(cellchat.second.art)
cellchat.second.art <- identifyOverExpressedInteractions(cellchat.second.art)

# project gene expression data onto protein-protein interaction (PPI)
cellchat.second.art <- projectData(cellchat.second.art, PPI.human) # PPI.mouse for mouse samples

# Compute the communication probability and infer cellular communication network
cellchat.second.art <- computeCommunProb(cellchat.second.art, raw.use = FALSE) 

# Filter out the cell-cell communication min.cells = 10
cellchat.second.art <- filterCommunication(cellchat.second.art, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat.second.art <- computeCommunProbPathway(cellchat.second.art)

# Calculate the aggregated cell-cell communication network
cellchat.second.art <- aggregateNet(cellchat.second.art)

# Compute the network centrality scores
cellchat.second.art <- netAnalysis_computeCentrality(cellchat.second.art, slot.name = "netP")

#' Figure 5A ====================================================================================================

first_IA<- netVisual_chord_gene(cellchat.first.art, 
                          sources.use = 3, 
                          targets.use = c(0,2, 6, 7, 11, 17, 18, 21), 
                          lab.cex = 1.2,
                          legend.pos.y = 5,
                     small.gap = 2, 
                     pairLR.use = data.frame('interaction_name'= c("MIF_CD74_CD44", "CCL5_CCR1", "CCL5_CCR5", 
                     "CCL5_CCR3", "CXCL9_CXCR3", "CXCL10_CXCR3", "CXCL11_CXCR3")),
                     show.legend = FALSE)
first_IA

Second_IA<- netVisual_chord_gene(cellchat.second.art,
                          sources.use = 3, 
                          targets.use = c(0,2, 6, 7, 11, 17, 18, 21), 
                          lab.cex = 1.2,
                          legend.pos.y = 5,
                     small.gap = 2, 
                     pairLR.use = data.frame('interaction_name'= c("MIF_CD74_CD44", "CCL5_CCR1", "CCL5_CCR5", 
                     "CCL5_CCR3", "CXCL9_CXCR3", "CXCL10_CXCR3", "CXCL11_CXCR3")),
                     show.legend = FALSE)


#' Figure 5B =====================================================================================================

first_IA<- netVisual_chord_gene(cellchat.first.art, 
                          sources.use = 8 
                          targets.use = c(0,2, 6, 7, 11, 17, 18, 21), 
                          lab.cex = 1.2,
                          legend.pos.y = 5,
                     small.gap = 2, 
                     pairLR.use = data.frame('interaction_name'= c("MIF_CD74_CD44", "CCL5_CCR1", "CCL5_CCR5", 
                     "CCL5_CCR3", "CXCL9_CXCR3", "CXCL10_CXCR3", "CXCL11_CXCR3")),
                     show.legend = FALSE)
first_IA

Second_IA<- netVisual_chord_gene(cellchat.second.art,
                          sources.use = 8 
                          targets.use = c(0,2, 6, 7, 11, 17, 18, 21), 
                          lab.cex = 1.2,
                          legend.pos.y = 5,
                     small.gap = 2, 
                     pairLR.use = data.frame('interaction_name'= c("MIF_CD74_CD44", "CCL5_CCR1", "CCL5_CCR5",     
                     "CCL5_CCR3", "CXCL9_CXCR3", "CXCL10_CXCR3", "CXCL11_CXCR3")),
                     show.legend = FALSE)



#' Figure 5C ======================================================================================================

first_IA_CXCL13<- netVisual_chord_gene(cellchat.first.art, 
                          sources.use = 10, 
                          targets.use = c(0, 3), 
                          lab.cex = 0.5,
                          legend.pos.y = 5,
                     small.gap = 2, 
                     signaling = c("CXCL", "CCL", "MIF", "TNF", "IL2", "IL1", "IL21"), 
                     show.legend = TRUE)
first_IA_CXCL13

second_IA_CXCL13<- netVisual_chord_gene(cellchat.second.art, 
                          sources.use = 10, 
                          targets.use = c(0, 3), 
                          lab.cex = 0.5,
                          legend.pos.y = 5,
                     small.gap = 2, 
                     signaling = c("CXCL", "CCL", "MIF", "TNF", "IL2", "IL1", "IL21"), 
                     show.legend = TRUE)
second_IA_CXCL13

# end of figure 5 ===========================================================================================






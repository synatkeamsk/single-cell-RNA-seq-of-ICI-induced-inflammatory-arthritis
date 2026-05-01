
library(Seurat)
library(patchwork)
library(sctransform)
library(glmGamPoi)
suppressMessages(require(tidyverse))
suppressMessages(library(scRepertoire))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

#' figure 2A ========================================================================================
sub.T.obj<- read_rds("umap.obj.pc15.rds")
DimPlot(sub.T.obj, 
reduction= "umap", 
label = TRUE, 
sizes.highlight = 10) +
  theme_minimal() +
  NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#' figure 2B ======================================================================================== 
overlay<- clonalOverlay(integrated.obj.T, 
              reduction = "umap", 
              bin= 20, 
              facet.by = "type") + 
  theme(plot.title = element_blank(), 
        legend.text = element_text(size= 13, face = "bold"), 
        strip.text = element_text(size = 14, face = "bold")) +
  NoAxes() + 
  NoLegend()
overlay
ggsave("umap.integrated.pdf", plot = overlay, height=5.5, width =10, units = "in", dpi = 300) 

#' figure 2C ======================================================================================== 
relativeabundance<- clonalOccupy(integrated.obj.T, 
                     x.axis = "ident", 
                     proportion = TRUE, 
                     label = FALSE, 
                     palette = "OrRd") + 
  theme_classic() + 
  theme(axis.text.x = element_text(hjust=0.5,vjust=0.2, face = "bold", size = 18), 
         axis.text.y = element_text(face = "bold", size = 18), 
         axis.title = element_text(face = "bold", size = 15),
         legend.title = element_blank(), 
         legend.text = element_text(size= 13, face = "bold"),
         legend.position = "top") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  labs(x= "Clusters") 
relativeabundance
ggsave("relative.abundance.pdf", plot = relativeabundance, height=6, width =8.5, units = "in", dpi = 300) 

#' figure 2D ======================================================================================== 
shareclone<- clonalOverlap(integrated.obj.T, 
              cloneCall = "strict",
              method = "morisita",
              group.by = "ident", 
              palette = "OrRd") + 
  theme(axis.text = element_text(size= 16, face = "bold"), 
        axis.title= element_text(size= 14, face = "bold")) + 
  labs(x= "Clusters", y= "Clusters")
shareclone
ggsave("heatmap.overalap.pdf", plot = shareclone, height=5, width =8.5, units = "in", dpi = 300) 

#' figure 2E ======================================================================================== 
alluvial_C164<- clonalCompare(C_9, 
                  top.clones = 8, 
                  group.by = "orig.ident",
                  samples = c("164S", "164S2"),   
                  cloneCall = "CTaa", 
                  graph = "alluvial", 
                  palette = "Dark 3", 
                  relabel.clones = TRUE) + 
  theme_classic() +
   theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
     axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2, face = "bold", size = 18), 
         axis.text.y = element_text(face = "bold", size = 18), 
         axis.title.y = element_text(face = "bold", size = 20), 
          axis.title.x = element_blank(),
         legend.title = element_text(size = 18, face = "bold"), 
         legend.text = element_text(size = 14, face = "bold"))
alluvial_C164
ggsave("alluvial_164.pdf", plot = alluvial_C164, height=5, width =5, units = "in", dpi = 300)

#' figure 2F ======================================================================================== 
  alluvial_C184<- clonalCompare(C_9, 
                  top.clones = 7, 
                  group.by = "orig.ident",
                  samples = c("184S2", "184S3"),   
                  cloneCall = "CTaa", 
                  graph = "alluvial", 
                  palette = "Dark 3", 
                  relabel.clones = TRUE) + 
   theme_classic() +
   theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
     axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2, face = "bold", size = 18), 
         axis.text.y = element_text(face = "bold", size = 18), 
         axis.title.y = element_text(face = "bold", size = 20), 
          axis.title.x = element_blank(),
         legend.title = element_text(size = 18, face = "bold"), 
         legend.text = element_text(size = 14, face = "bold"))
alluvial_C184
ggsave("alluvial_184.pdf", plot = alluvial_C184, height=5, width =5.5, units = "in", dpi = 300) 

#' figure 2G ======================================================================================== 
alluvial_218<- clonalCompare(C_9, 
                  top.clones = 10, 
                  group.by = "orig.ident",
                  samples = c("218S", "218S2"),   
                  cloneCall = "CTaa", 
                  graph = "alluvial", 
                  relabel.clones = TRUE, 
                  palette = "Dark 3") + 
theme_classic() +
   theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
     axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2, face = "bold", size = 18), 
         axis.text.y = element_text(face = "bold", size = 18), 
         axis.title.y = element_text(face = "bold", size = 20), 
          axis.title.x = element_blank(),
         legend.title = element_text(size = 18, face = "bold"), 
         legend.text = element_text(size = 14, face = "bold"))
alluvial_218
ggsave("alluvial_218.pdf", plot = alluvial_218, height=5, width =5, units = "in", dpi = 300) 

#' figure 2H ========================================================================================
alluvial_164C3<- clonalCompare(C_3, 
                  top.clones = 10, 
                  group.by = "orig.ident",
                  samples = c("164S", "164S2"),   
                  cloneCall = "CTaa", 
                  graph = "alluvial", 
                  relabel.clones = TRUE, 
                  palette = "Dark 3") +
theme_classic() +
   theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
     axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2, face = "bold", size = 18), 
         axis.text.y = element_text(face = "bold", size = 18), 
         axis.title.y = element_text(face = "bold", size = 20), 
          axis.title.x = element_blank(),
         legend.title = element_text(size = 18, face = "bold"), 
         legend.text = element_text(size = 14, face = "bold"))
alluvial_164C3
ggsave("alluvial_164C3.pdf", plot = alluvial_164C3, height=5, width =5, units = "in", dpi = 300)

#' figure 2I ========================================================================================
alluvial_184C3<- clonalCompare(C_3, 
                  top.clones = 10, 
                  group.by = "orig.ident",
                  samples = c("184S2", "184S3"),   
                  cloneCall = "CTaa", 
                  graph = "alluvial", 
                  relabel.clones = TRUE, 
                  palette = "Dark 3")  +
  theme_classic() +
   theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
     axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2, face = "bold", size = 18), 
         axis.text.y = element_text(face = "bold", size = 18), 
         axis.title.y = element_text(face = "bold", size = 20), 
          axis.title.x = element_blank(),
         legend.title = element_text(size = 18, face = "bold"), 
         legend.text = element_text(size = 14, face = "bold"))
alluvial_184C3
ggsave("alluvial_184C3.pdf", plot = alluvial_184C3, height=5, width =5, units = "in", dpi = 300) 

#' figure 2J ========================================================================================
alluvial_218C3<- clonalCompare(C_3, 
                  top.clones = 10, 
                  group.by = "orig.ident",
                  samples = c("218S", "218S2"),   
                  cloneCall = "CTaa", 
                  graph = "alluvial", 
                  relabel.clones = TRUE, 
                  palette = "Dark 3")  +
  theme_classic() +
   theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
     axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2, face = "bold", size = 18), 
         axis.text.y = element_text(face = "bold", size = 18), 
         axis.title.y = element_text(face = "bold", size = 20), 
          axis.title.x = element_blank(),
         legend.title = element_text(size = 18, face = "bold"), 
         legend.text = element_text(size = 14, face = "bold"))
alluvial_218C3
ggsave("alluvial_218C3.pdf", plot = alluvial_218C3, height= 5, width =5, units = "in", dpi = 300)

#' end of figure 2 ===============================================================================================

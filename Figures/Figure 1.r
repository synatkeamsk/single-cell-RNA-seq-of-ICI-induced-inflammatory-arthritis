library(Seurat)
library(tidyverse)
library(patchwork)
library(viridis)

#Load Seurat preprocessed obj
seurat.obj<- readRDS("seurat.obj.final.RDS")

#'  figure 1C ========================================================================================================
DimPlot(seurat.obj, reduction = "umap", label = FALSE) +
  theme_classic() +
  NoLegend() + 
  theme(axis.title = element_text(size= 15, face = "bold"),
        axis.text = element_text(size = 14, face = "bold")) + 
  NoAxes()

#' figure 1D =========================================================================================================
genes <- c("CD3D", "MKI67", "CD4", "IL2RA", "FOXP3", "CD8A", "TRDV2", "IL7R", "KLRG1", "NCAM1", "CD19", "TBX21", "MS4A1", "CD27", "CD38", "SDC1", "CD14", 
"FCGR3A", "S100A8", "S100A9", "HLA-DRB1", "CD1C", "CXCR3", "CXCR6", "CCR6", "CCR7", "CXCL13", "TOX", "LAG3", "PDCD1","CTLA4", "HAVCR2", "TGIT", "EOMES", 
"GZMB", "GZMH", "GZMK", "IFNG", "TNF", "IL10", "IL17A", "IL21", "TGFB1", "IL2")

# Define the order of clusters (from top to bottom)
desired_order <- c("CD4 Tcm", "Recently activated CD4", "PD1hi CXCL13hi CD4", "Treg", "Effector CD8", 
                   "CD8 Trm", "Cycing T", "MAIT", "Mixed T", "NK", "NK T", "mDC", "CLEC9A+ DCs", 
                   " pDCs", "Non-inflammatory DCs", "Classical mono", "Non-classical mono", "SPP1+ Mac", 
                   "Neutrophil", "B cells", "Syn cells")

# Reorder the identities in Seurat object
all.annotated.subsets <- SetIdent(all.annotated.subsets, 
                                  value = factor(Idents(all.annotated.subsets), 
                                                 levels = desired_order))

Dot.plot<- DotPlot(all.annotated.subsets, features = genes, dot.scale = 4.5) + 
 geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="inferno") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15), 
        axis.title.y= element_blank(), 
        axis.text = element_text(size=10, face = "bold"), 
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2, face = "bold", size= 10), 
        axis.title.x = element_blank(),
        legend.title = element_text(size= 12, face = "bold"), 
        legend.text = element_text(size= 10, face = "bold"), 
        legend.position = "right") 

#' figure 1E  =========================================================================================================
proportion_normalized <- proportion %>%
  group_by(Group) %>%
  mutate(percent = frequency / sum(frequency) * 100)

# order group
proportion_normalized$Group<- factor(proportion_normalized$Group, levels = 
                                       c("osteoarthritis", "ICI_arthritis"))
#proportion_normalized <- subset(proportion_normalized, !is.na(cluster_name))


proportion_normalized$cluster_name<- factor(proportion_normalized$cluster_name, levels = 
                                       c("Neutrophil", 
                                         "Classical Monocytes", 
                                         "Non-classical monocytes", 
                                         "SPP1+ Macrophages", 
                                         "mDC", 
                                         "pDC", 
                                         "CD1c+ Non-inflammatory DCs", 
                                         "CLEC9A+ DCs", 
                                         "CD4 central memory", 
                                         "recently activated CD4", 
                                         "PD-1hi CXCL13hi CD4", 
                                         "Treg", 
                                        "CD8 Trm", 
                                         "Effector CD8", 
                                         "Cycing T", 
                                         "MAIT", 
                                         "Mixed T", 
                                         "NK", 
                                         "NKT",
                                         "B cells", 
                                         "Synovial cells"))

# ggplot
allimmune.stack<- ggplot(proportion_normalized, aes(x = Group, y = percent, fill = cluster_name)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  theme(axis.title = element_blank(), 
        axis.text.y = element_text(size=8, face= "bold"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        legend.text = element_text(size = 8), 
        legend.title = element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
  scale_fill_manual(values = c(
    "Effector CD8" = "#FF0000", 
    "PD-1hi CXCL13hi CD4"= "#0000FF", 
    "CD4 central memory" = "#30D5C8", 
    "Treg"= "#FFC0CB", 
    "Non-classical monocytes"= "#F8766D", 
    "SPP1+ Macrophages"= "#ED8141",
    "mDC"= "#CF9400", 
    "pDC"= "#BB9D00", 
    "CD1c+ Non-inflammatory DCs"= "#95A", 
    "CLEC9A+ DCs"="#00B81F", 
    "recently activated CD4"= "#00BF7D",
    "Cycing T"= "#00ABFD", 
    "MAIT"= "#9590FF", 
    "Mixed T"= "#DC71FA", 
    "NK"= "#FF62BC", 
    "NKT"= "#00C1AA",
    "B cells"= "#FC717F", 
    "Synovial cells"= "#000000",
    "Neutrophil"= "#AC88FF", 
    "Classical Monocytes"= "#00BE6C", 
    "CD8 Trm" ="#808080")) 
allimmune.stack

#' figure 1E =========================================================================================================
CD4.ctm<- proportion %>% filter(cluster_name == "CD4 central memory")
CD4.ctm$Group<- factor(CD4.ctm$Group, levels = c("osteoarthritis", "ICI_arthritis"))
CD4.ctm$Arthritis<- factor(CD4.ctm$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))

#visualization ! 
p1<- ggplot(CD4.ctm, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(), 
  legend.position = "none") + 
  labs(x= "Groups", y= "% of Total Cells", title = "CD4 Tcm")
p1

#PD1hiCXCL13hi CD4
TPH<- proportion %>% filter(cluster_name == "PD-1hi CXCL13hi CD4")
TPH$Group<- factor(TPH$Group, levels = c("osteoarthritis", "ICI_arthritis"))
TPH$Arthritis<- factor(TPH$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))


#visualization ! 
p3<- ggplot(TPH, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face = "bold"), 
  axis.text.x = element_blank(),
  legend.position = "none") + 
  labs(x= "Groups", y= "% of Total Cells", title = "PD1hi CXCL13hi CD4")
p3

#Regulatory T
Treg<- proportion %>% filter(cluster_name == "Treg")
Treg$Group<- factor(Treg$Group, levels = c("osteoarthritis", "ICI_arthritis"))
Treg$Arthritis<- factor(Treg$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))


#visualization ! 
p5<- ggplot(Treg, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face= "bold"), 
  axis.text.x = element_blank(),
  legend.position = "none") + 
  labs(x= "Groups", y= "% Treg of total cells", title = "Treg")
p5

# Effector CD8
effector.cd8<- proportion %>% filter(cluster_name == "Effector CD8")
effector.cd8$Group<- factor(effector.cd8$Group, levels = c("osteoarthritis", "ICI_arthritis"))
effector.cd8$Arthritis<- factor(effector.cd8$Arthritis, levels = c("osteoarthritis", "first arthritis", "second arthritis"))


#visualization ! 
p6<- ggplot(effector.cd8, aes(x= Group, y= frequency, fill= Group)) +
  theme_classic() +
  geom_boxplot(size= 0.5, alpha= 0.5) +
   geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=Arthritis), show.legend = TRUE, size= 2.5) +
 theme(plot.title = element_text(hjust = 0.5, face = "bold", size=5.5
  ), axis.title.y = element_blank(), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=8, face = "bold"), 
  axis.text.x = element_blank(),
  legend.position = "none") + 
  labs(x= "Groups", y= "% of Total Cells", title = "Effector CD8")
p6

combine<- (p1|p3|p5|p6)
ggsave("immune_fraction.rna.pdf", plot = combine, height=1.5, width = 4.5, units = "in", dpi = 300)

#' figure 1F ============================================================================================
flow.data<- read.csv("Tcell.csv", 
                   stringsAsFactors = TRUE, 
                   fileEncoding = "ISO-8859-1", 
                   check.names = FALSE)

PD1hiCXCL13hiCD4<- ggplot(flow.data, aes(x = Group, y = PD1hiCXCL13hiCD4, 
                                       fill= Group)) +
  theme_classic() +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=type), show.legend = TRUE, size= 1.7) +
  scale_shape_manual(values = c(16, 17, 18, 19)) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.text.x = element_blank()) + 
  theme(legend.position = "none") + 
  labs(x = "Groups", y = "% of live cells", title = "PD1hi CXCL13hi CD4+ T cells") 
PD1hiCXCL13hiCD4


CD4Tcm<- ggplot(heatmap, aes(x = Group, 
                                   y = CD4Tcm, 
                             fill= Group
)) +
  theme_classic() +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=type), show.legend = TRUE, size= 1.7) +
  scale_shape_manual(values = c(16, 17, 18, 19)) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.text.x = element_blank()) + 
  theme(legend.position = "none") + 
  labs(x = "Groups", y = "% of live cells", title = "CD4 Central memory") 
CD4Tcm

Treg<- ggplot(heatmap, aes(x = Group, y = Treg, 
                           fill= Group)) +
  theme_classic() +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=type), show.legend = TRUE, size= 1.7) +
  scale_shape_manual(values = c(16, 17, 18, 19)) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.text.x = element_blank()) + 
  theme(legend.position = "none") + 
  labs(x = "Groups", y = "% of live cells", title = "Tregs") 
Treg

effectorcd8<- ggplot(heatmap, aes(x = Group, y = EffectorCD8, fill= Group)) +
  theme_classic() +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(color=type), show.legend = TRUE, size= 1.7) +
  scale_shape_manual(values = c(16, 17, 18, 19)) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.text.x = element_blank()) + 
  theme(legend.position = "none") + 
  labs(x = "Groups", y = "% of live cells", title = "Effector CD8") 
effectorcd8

combine<- (CD4Tcm|PD1hiCXCL13hiCD4|Treg|effectorcd8)
ggsave("immune_fraction.annals.flow.pdf", plot = combine, height=1.5, width = 4.5, units = "in", dpi = 300)

#' end of figure 1 ==========================================================================================
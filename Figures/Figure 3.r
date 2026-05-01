
library(tidyverse) 
library(RColorBrewer) 
library(ggrepel) 
library(Harmony)
library(monocle3)
library(SeuratWrappers)

seurat.obj.T<- readRDS("umap.obj.pc15.rds")
DimPlot(seurat.obj.T, label = TRUE) + 
  NoAxes() + 
  NoLegend() 

T.annotated.subsets<- RenameIdents(seurat.obj.T, `0` ="Th17 like",`10`= "Activated CD4", `4`= "CD4 Tcm", 
                                   `6`= "PD-1hi CXCL13hi CD4", `8`= "Th1-like", `15`= "MAIT", `3` ="Treg", `11`="Treg", 
                                   `18`= "MAIT", `2` ="gd T", `7`= "Mixed T", `5`= "CD8 Tem", `1` ="CD8 effector",
                                   `14`= "CD8 Trm",  `9`= "CD8 effector",  `17`= "KLRB1hi CD8 Tem",  `12`= "NKG7+ CD8", `16`= "Cycling T", 
                                   `13`= "NKT")
T.annotated.subsets@meta.data$group <- ifelse(T.annotated.subsets@meta.data$type %in% c("First arthritis", "Second arthritis"), "ICI-arthritis", "Osteoarthritis")

# Identify differentially expressed genes across conditions
T.annotated.subsets$celltypes<-Idents(T.annotated.subsets)
T.annotated.subsets$arthritis_celltype <- paste(T.annotated.subsets$group, sep = "_", Idents(T.annotated.subsets))
Idents(T.annotated.subsets) <- "arthritis_celltype"
prepsct<- PrepSCTFindMarkers(T.annotated.subsets)
effector_cd8_icivsosteo<- FindMarkers(T.annotated.subsets,assay = "SCT", ident.1 = "ICI-arthritis_CD8 effector", 
                               ident.2 = "Osteoarthritis_CD8 effector", min.pct = 0.25, logfc.threshold = 0.25, 
                           verbose = FALSE, 
                           recorrect_umi = FALSE)

head(effector_cd8_icivsosteo, n = 20)

# order it from highest log2FC to lowest
effector_cd8_icivsosteo_order<- effector_cd8_icivsosteo %>% arrange(desc(avg_log2FC))
head(effector_cd8_icivsosteo_order)
effector_cd8_icivsosteo_order<-effector_cd8_icivsosteo_order[effector_cd8_icivsosteo_order$p_val_adj<0.05, ]
# Save deg of cluster 1
write.csv(effector_cd8_icivsosteo_order, file = "degicivsosteo.effector.CD8.csv")

#read back the file
deg.ici.vs.osteo<- read.csv("deg_between conditions/deg_ici_vs_osteo.effector.CD8.csv", row.names = 1)
head(deg.ici.vs.osteo)

#' figure 3A =======================================================================================================
volcano_effectorcd8icivsoste<- ggplot(data = deg.ici.vs.osteo,
                                      aes(x = avg_log2FC, y = -log10(p_val_adj), 
                                                                   col = diffexpressed, 
                                          label= delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  theme(plot.title = element_text(hjust = 0.5, size = 8.5, face = "bold"), 
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(face = "bold", size=11), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 7.3), 
        legend.position = c(0.85, 0.23)) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), 
                     labels = c("Downregulated", "Not significant", "Upregulated"))  + 
 ggtitle("ICI-arthritis vs. Osteoarthritis (effector CD8 cluster)") + 
   geom_text_repel(max.overlaps = Inf, size=2.3) # To show all labels 
volcano_effectorcd8icivsoste

# save high quality plot! 
ggsave("volcanoiciosteo_effectorcd8.new.pdf", plot = volcano_effectorcd8icivsoste, height=4, width = 4, units = "in", dpi = 300)

# Gene Set Enrichment analysis 
effector_cd8_icivsosteo$gene<- rownames(effector_cd8_icivsosteo)
effector_cd8_icivsosteo<- effector_cd8_icivsosteo %>% arrange(desc(avg_log2FC))
fold_changes<- effector_cd8_icivsosteo$avg_log2FC
names(fold_changes)<- effector_cd8_icivsosteo$gene

#Load geneset ! 
Reactome <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
hallmark <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/h.all.v2023.1.Hs.symbols.gmt")
KEGG <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
GO <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c5.go.bp.v2023.1.Hs.symbols.gmt")

#combine the four gsea
Four_gene_sets <- c(hallmark, KEGG, GO, Reactome)

#getset for all pathway 
Allgenesets <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/msigdb.v2023.1.Hs.symbols.gmt")

# GSEA analysis
gsea<- fgsea(pathways = hallmark,
                  stats = fold_changes,
                  eps = 0.0,
                  minSize=15,
                  maxSize=500)

head(gsea_Myo[order(pval), ]) 

#save csv file! 
write.csv(gsea_Myo, file = "gsea_ici_osteo.csv")

#' figure 3B =======================================================================================================
gseaici_oste<- read.csv("gseaici_osteo.csv", stringsAsFactors = TRUE, header = TRUE)
gseaici<- ggplot(gseaici_oste, aes(reorder(Pathway, NES), NES, fill= p.adjust)) + 
  theme_classic() +
  geom_bar(stat = "identity")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y  = element_blank(), 
        axis.text = element_text(face = "bold", size = 10), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10)) + 
  labs(x="Hallmark pathways", y="normalized enrichment score",
       title="ICI-arthritis vs. Osteoarthritis (Effector CD8)") + 
  coord_flip() + 
  scale_fill_gradient(low = "red", high = "blue") + 
   NoAxes() + 
  NoLegend()
gseaici
ggsave("gseaici.pdf", plot = gseaici, height=3, width =6.5, units = "in", dpi = 300)

#' figure 3C =======================================================================================================
T.annotated.subsets$celltypes<-Idents(T.annotated.subsets)
T.annotated.subsets$arthritis_celltype <- paste(T.annotated.subsets$type, sep = "_", Idents(T.annotated.subsets))
Idents(T.annotated.subsets) <- "arthritis_celltype"

#Calculate differentially expressed genes between first and second arthritis!
prepsct<- PrepSCTFindMarkers(T.annotated.subsets)
effector_cd8<- FindMarkers(T.annotated.subsets,assay = "SCT", ident.1 = "Second arthritis_CD8 effector", 
                               ident.2 = "First arthritis_CD8 effector", min.pct = 0.25, logfc.threshold = 0.25, 
                           verbose = FALSE, 
                           recorrect_umi = FALSE)

# order it from highest log2FC to lowest
effector_cd8_order<- effector_cd8 %>% arrange(desc(avg_log2FC))
effector_cd8_order<- effector_cd8_order[effector_cd8_order$p_val_adj<0.05, ]
write.csv(effector_cd8_order, file = "deg.effector.CD8.csv")

#read the deg data
deg.first.vs.second.effectorCD8<- read.csv("deg_between conditions/deg_second_vs_first_effector.CD8.csv", row.names = 1)
deg.first.vs.second.effectorCD8$diffexpressed<- "NO"
deg.first.vs.second.effectorCD8$diffexpressed[deg.first.vs.second.effectorCD8$avg_log2FC > 0.5 & deg.first.vs.second.effectorCD8$p_val_adj < 0.05] <- "UP"
deg.first.vs.second.effectorCD8$diffexpressed[deg.first.vs.second.effectorCD8$avg_log2FC < -0.5 & deg.first.vs.second.effectorCD8$p_val_adj < 0.05] <- "DOWN"
head(deg.ici.vs.osteo[order(deg.first.vs.second.effectorCD8$p_val_adj) & deg.first.vs.second.effectorCD8$diffexpressed == 'DOWN', ])

#Volcano plot
volcano_effectorcd8_secondvsfirst<- ggplot(data = deg.first.vs.second.effectorCD8,
                                      aes(x = avg_log2FC, y = -log10(p_val_adj), 
                                                                   col = diffexpressed, 
                                          label= delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  theme(plot.title = element_text(hjust = 0.5, size = 8.5, face = "bold"), 
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(face = "bold", size=11), 
        legend.title = element_blank(), 
        legend.text = element_blank(), 
        legend.position = ) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), 
                     labels = c("Downregulated", "Not significant", "Upregulated"))  + 
 ggtitle("Second vs. first arthritis (effector CD8 cluster)") + 
   geom_text_repel(max.overlaps = Inf, size=2.3) # To show all labels 
volcano_effectorcd8_secondvsfirst
ggsave("volcano_secondvsfirst_effectorcd8.pdf", plot = volcano_effectorcd8_secondvsfirst, height=4, width = 4.9, units = "in", dpi = 300)


#' figure 3D =======================================================================================================
effector_cd8$gene<- rownames(effector_cd8)
effector_cd8<- effector_cd8 %>% arrange(desc(avg_log2FC))
fold_changes<- effector_cd8$avg_log2FC
names(fold_changes)<- effector_cd8$gene
Reactome <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
hallmark <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/h.all.v2023.1.Hs.symbols.gmt")
KEGG <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
GO <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c5.go.bp.v2023.1.Hs.symbols.gmt")
Four_gene_sets <- c(hallmark, KEGG, GO, Reactome)

#getset for all pathway 
Allgenesets <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/msigdb.v2023.1.Hs.symbols.gmt")
gsea <- fgsea(pathways = hallmark,
                  stats = fold_changes,
                  eps = 0.0,
                  minSize=15,
                  maxSize=500)

#read back the csv file of GSEA
gseaicisecondvsfirst<- read.csv("enrichment.csv", stringsAsFactors = TRUE, header = TRUE)

# ggplot
gseaicisecondvsfirst<- ggplot(gseaicisecondvsfirst, aes(reorder(Pathway, NES), NES, fill= padj)) + 
  theme_classic() +
  geom_bar(stat = "identity")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y  = element_blank(), 
        axis.text = element_text(face = "bold", size = 10), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10)) + 
  labs(x="Hallmark pathways", y="normalized enrichment score",
       title="Second vs. first arthritis (Effector CD8)") + 
  coord_flip() + 
  scale_fill_gradient(low = "red", high = "blue")
gseaicisecondvsfirst
ggsave("gseasecondvsfirst.pdf", plot = gseaicisecondvsfirst, height=3, width =6, units = "in", dpi = 300)

#' figure 3F =======================================================================================================
cytokine<- read.csv("cytokine_T.csv", 
                    fileEncoding = "ISO-8859-1")
#order the factor !! 
cytokine$type<- factor(cytokine$type, levels = c("OA", "ICI-IA"))
cytokine$group<- factor(cytokine$group, levels = c("first RA", "second RA"))

#IFy+ TNFa-
Cd8ifngpos_tnfneg<- ggplot(firstvssecond, aes(group,TNFnegIFNgpos, 
                                              shape = group, fill= group)) +
  geom_boxplot(alpha=0.5) +
  theme_classic() +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(shape = group), show.legend = TRUE, size= 1.8) +
  scale_shape_manual(values = c(16, 17)) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_blank()) +
  theme(legend.position = "none")  + 
  labs(x = "Groups", y = "% of CD8", title = "TNF-A")
Cd8ifngpos_tnfneg
ggsave("cd8ifng1.pdf", plot = Cd8ifngpos_tnfneg, height=1.2, width = 1.5, units = "in", dpi = 300)

# TNFpos IFNG- 
Cd8tnfpos_ifngneg<- ggplot(firstvssecond, aes(group,TNFposIFNgneg, 
                                   shape = group, fill= group)) +
  geom_boxplot(alpha=0.5) +
  theme_classic() +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(shape = group), show.legend = TRUE, size= 1.8) +
  scale_shape_manual(values = c(16, 17)) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_blank()) +
  theme(legend.position = "none")  + 
  labs(x = "Groups", y = "% of CD8", title = "TNF-A")
Cd8tnfpos_ifngneg
ggsave("cd8tnfa.pdf", plot = Cd8tnfpos_ifngneg, height=1.2, width = 1.5, units = "in", dpi = 300)

#TNFa+IFNy+ 
effectorcd8tnfifng<- ggplot(firstvssecond, aes(group,TNFaIFNg_of_CD8, fill=group, 
                                               shape = group)) +
  geom_boxplot(alpha=0.5) +
  theme_classic() +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point( aes(fill=group,group=paired), position = position_dodge(0.2), 
             size= 1.8) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_blank()) +
  theme(legend.position = "none")  + 
  labs(x = "Groups", y = "% of CD8", title = "TNF / IFNG CD8")
effectorcd8tnfifng
ggsave("cd8ifngtnfa.pdf", plot = effectorcd8tnfifng, height=1.2, width = 1.5, units = "in", dpi = 300)

#' figure 3H =======================================================================================================
TC1<- ggplot(firstvssecond, aes(group,Tc1_of_CD8, 
                                             shape = group, fill= group)) +
  geom_boxplot(alpha=0.5) +
  theme_classic() +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(shape = group), show.legend = TRUE, size= 3) +
  scale_shape_manual(values = c(16, 17)) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_blank()) +
  theme(legend.position = "none")  + 
  labs(x = "Groups", y = "% of CD4", title = "TC1")
TC1
ggsave("TC1.pdf", plot = TC1, height=2, width = 2, units = "in", dpi = 300)

TC17<- ggplot(firstvssecond, aes(group,Tc17_of_CD8, 
                                shape = group, fill= group)) +
  geom_boxplot(alpha=0.5) +
  theme_classic() +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(shape = group), show.legend = TRUE, size= 3) +
  scale_shape_manual(values = c(16, 17)) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_blank()) +
  theme(legend.position = "none")  + 
  labs(x = "Groups", y = "% of CD4", title = "TC1")
TC17
ggsave("TC17.pdf", plot = TC17, height=2, width = 2, units = "in", dpi = 300)

#Tc1/17 
TC1.17<- ggplot(firstvssecond, aes(group,Tc1.17_of_CD8, 
                                 shape = group, fill= group)) +
  geom_boxplot(alpha=0.5) +
  theme_classic() +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(position=position_jitterdodge(jitter.width=0.3, dodge.width = 0.3), 
             aes(shape = group), show.legend = TRUE, size= 3) +
  scale_shape_manual(values = c(16, 17)) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_blank()) +
  theme(legend.position = "none")  + 
  labs(x = "Groups", y = "% of CD4", title = "TC1.17") + 
  ylim(c(0,2.5))
TC1.17
ggsave("TC1.17.pdf", plot = TC1.17, height=2, width = 2, units = "in", dpi = 300)

#' figure 3I =======================================================================================================
seurat.obj.CD8<- subset(sub.T.obj, idents = c("1", "5", "9", 
                                                    "12", "14"))
DimPlot(seurat.obj.CD8, reduction= "umap", label = TRUE) +
  theme_minimal() +
  NoLegend()

# Run SCT/PCA 
cd8.subcluster<- SCTransform(seurat.obj.CD8, 
                         method= "glmGamPoi")
cd8.subcluster<- RunPCA(cd8.subcluster, verbose = FALSE)
ElbowPlot(cd8.subcluster, ndims = 50,reduction = "pca") +
  ggtitle("Number of principle component")  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
cd8.T.integration<- RunHarmony(cd8.subcluster,
                              group.by.vars = "orig.ident")
CD8.T.pc15 <- FindNeighbors(cd8.T.integration, 
                                 reduction = "harmony", 
                               dims = 1:15)
CD8.T.pc15<- FindClusters(CD8.T.pc15,
                                 resolution = 0.11)
CD8.T.pc15<- RunUMAP(object = CD8.T.pc15,
                          reduction = "harmony",
                          dims = 1:15)
cd8_subclu= DimPlot(CD8.T.pc15, 
                    reduction= "umap", 
                    label = FALSE, 
                    sizes.highlight = 10) +
  theme_classic() +
  NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        panel.grid = element_blank(),
        axis.text = element_text(face = "bold", size= 16), 
        axis.title = element_text(face = "bold", size= 18)) + 
  NoAxes()
cd8_subclu
ggsave("CD8.new.pdf", plot = cd8_subclu, height=4, width = 4.5, units = "in", dpi = 300)

#' figure 3J =======================================================================================================
seurat_obj_cd8<- readRDS("cd8_subclu.rds")

#Dimplot 
DimPlot(seurat_obj_cd8, reduction = "umap", label = TRUE) +
  theme_classic() +
  NoLegend() + 
  NoAxes()

# Second arthritis 
second_arthritis<- subset(seurat_obj_cd8, subset = type == c("Second arthritis"))
monocle.cds.s <- as.cell_data_set(second_arthritis)
monocle.cds.s <- cluster_cells(cds = monocle.cds.s, reduction_method = "UMAP")
monocle.cds.s <- learn_graph(monocle.cds.s)
fData(monocle.cds.s) 
rownames(fData(monocle.cds.s))[1:10]
fData(monocle.cds.s)$gene_short_name <- rownames(fData(monocle.cds.s))
fData(monocle.cds.s)
startCells <- Cells(second_arthritis)[second_arthritis$seurat_clusters == 3]
monocle.cds.s <- order_cells(monocle.cds.s, 
                           reduction_method = "UMAP",
                           root_cells = startCells)
cd8_pseudo<- plot_cells(cds = monocle.cds.s,
           color_cells_by = "pseudotime",
           show_trajectory_graph = T,
           trajectory_graph_color = "red",
           trajectory_graph_segment_size = 1,
           graph_label_size = 2,
           cell_size = 0.8,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F) + 
  theme_bw() + 
theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(face = "bold", size= 12), 
        axis.title = element_text(face = "bold", size= 12), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) + 
  theme(legend.position = "none")
cd8_pseudo

#plot the pseudotime !!
monocle.cds.s$monocle3_pseudotime <- pseudotime(monocle.cds.s)
data.pseudo_second <- as.data.frame(colData(monocle.cds.s))
pseudotime_second<- ggplot(data.pseudo_second, aes(x = monocle3_pseudotime, y = reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) +
  geom_density_ridges() + 
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title.y  = element_blank(), 
        axis.text = element_text(size= 12, face = "bold"), 
        axis.title.x = element_text(siz= 16, face = "bold"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  labs(x= "Pseudotime")
pseudotime_second
ggsave("cd8.pseudo.2.pdf", plot = pseudotime_second, height=3, width= 2.5,  units = "in", dpi = 300)
 panel.border = element_rect(colour = "black", fill = NA, size = 1.3)
ggsave("CD8secondtraje.pdf", plot = pseudotime_second, height=3, width = 3.5, units = "in", dpi = 300)

# First arthritis 
first_arthritis<- subset(seurat_obj_cd8, subset = type == "First arthritis")
monocle.cds.f <- as.cell_data_set(first_arthritis)
monocle.cds.f <- cluster_cells(cds = monocle.cds.f, 
                             reduction_method = "UMAP")
monocle.cds.f <- learn_graph(monocle.cds.f)
startCells <- Cells(first_arthritis)[first_arthritis$seurat_clusters == 3]
monocle.cds.f <- order_cells(monocle.cds.f,
                           reduction_method = "UMAP", 
                           root_cells = startCells)
plot_cells(cds = monocle.cds.f,
           color_cells_by = "pseudotime",
           show_trajectory_graph = T,
           trajectory_graph_color = "grey",
           trajectory_graph_segment_size = 0.5,
           graph_label_size = 2,
           cell_size = 1,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F)

#plot the pseudotime !!
monocle.cds.f$monocle3_pseudotime <- pseudotime(monocle.cds.f)
data.pseudo_first <- as.data.frame(colData(monocle.cds.f))

pseudotime_first<- ggplot(data.pseudo_first, aes(x = monocle3_pseudotime, y = reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) +
  geom_density_ridges() + 
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title.y  = element_blank(), 
        axis.text = element_text(size= 12, face = "bold"), 
        axis.title.x = element_text(siz= 16, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  labs(x= "Pseudotime")
pseudotime_first
ggsave("cd8.pseudo.1.pdf", plot = pseudotime_first, height=3, width= 3.5,  units = "in", dpi = 300)

#'   end-of-figure3 ========================================================================================
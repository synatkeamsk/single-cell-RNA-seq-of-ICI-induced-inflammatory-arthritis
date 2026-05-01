library(tidyverse) 
library(RColorBrewer) 
library(ggrepel) 
library(Harmony)
library(monocle3)
library(SeuratWrappers)
library(ggrepel)
library(tidyr)
library(ggplot2)
library(ggridges)

seurat.obj.T<- readRDS("umap.obj.pc15.rds")
DimPlot(seurat.obj.T, label = TRUE) + 
  NoAxes() + 
  NoLegend() 

T.annotated.subsets<- RenameIdents(seurat.obj.T, `0` ="Th17 like",`10`= "Activated CD4", `4`= "CD4 Tcm", 
                                   `6`= "PD-1hi CXCL13hi CD4", `8`= "Th1-like", `15`= "MAIT", `3` ="Treg", `11`="Treg", 
                                   `18`= "MAIT", `2` ="gd T", `7`= "Mixed T", `5`= "CD8 Tem", `1` ="CD8 effector",
                                   `14`= "CD8 Trm",  `9`= "CD8 effector",  `17`= "KLRB1hi CD8 Tem",  `12`= "NKG7+ CD8", `16`= "Cycling T", 
                                   `13`= "NKT")

T.annotated.subsets@meta.data$group <- ifelse(T.annotated.subsets@meta.data$type 
                                              %in% c("First arthritis", "Second arthritis"), "ICI-arthritis", "Osteoarthritis")

T.annotated.subsets$celltypes<-Idents(T.annotated.subsets)
T.annotated.subsets$arthritis_celltype <- paste(T.annotated.subsets$group, sep = "_", Idents(T.annotated.subsets))
Idents(T.annotated.subsets) <- "arthritis_celltype"

#Calculate differentially expressed genes between first and second arthritis!
prepsct<- PrepSCTFindMarkers(T.annotated.subsets)
pd1hicxcl13hicd4_icivsosteo<- FindMarkers(T.annotated.subsets,assay = "SCT", ident.1 = "ICI-arthritis_PD-1hi CXCL13hi CD4", 
                                      ident.2 = "Osteoarthritis_PD-1hi CXCL13hi CD4", min.pct = 0.25, logfc.threshold = 0.25, 
                                      verbose = FALSE, 
                                      recorrect_umi = FALSE)
head(pd1hicxcl13hicd4_icivsosteo)

# order it from highest log2FC to lowest
pd1hicxcl13hicd4_icivsosteo_order<- pd1hicxcl13hicd4_icivsosteo %>% arrange(desc(avg_log2FC))
head(pd1hicxcl13hicd4_icivsosteo_order)


pd1hicxcl13hicd4_icivsosteo_order<- pd1hicxcl13hicd4_icivsosteo_order[pd1hicxcl13hicd4_icivsosteo_order$p_val_adj<0.05, ]

# Save deg of cluster 1
write.csv(pd1hicxcl13hicd4_icivsosteo_order, file = "deg.cxcl13cd4icivs.osteo.csv")

volcano_pd1hicxcl13hicd4icivs.osteo <- pd1hicxcl13hicd4_icivsosteo[order(pd1hicxcl13hicd4_icivsosteo$p_val_adj),]

# Create data frame with cut_off based on significance and avg_log2FC thresholds
results_cxcl13icivs.osteo <- as.data.frame(mutate(as.data.frame(volcano_pd1hicxcl13hicd4icivs.osteo), 
                                    cut_off = case_when(
                                      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Downregulated",
                                      p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Upregulated",
                                      TRUE ~ "Not significant")), 
                             row.names = rownames(volcano_pd1hicxcl13hicd4icivs.osteo))

# Ensure the cut_off column is a factor with the correct levels
results_cxcl13icivs.osteo$cut_off <- factor(results_cxcl13icivs.osteo$cut_off, levels = c('Downregulated', 'Upregulated', 'Not significant'))

#' figure 4A =======================================================================================================
options(ggrepel.max.overlaps = 20)
volcano_cxcl13icivs.osteoplot <- ggplot(results_cxcl13icivs.osteo, aes(avg_log2FC, -log10(p_val_adj))) + 
  theme_classic() +
  geom_point(aes(col = cut_off), size= 2) +
  scale_color_manual(values = c("#00AFBB","#bb0c00", "grey")) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x= element_blank(),
        axis.text = element_text(face = "bold", size = 13), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 6.5), 
        legend.position = c(0.2, 0.85)) + 
  geom_text_repel(data = results_cxcl13icivs.osteo[1:30,], 
                  aes(label = rownames(results_cxcl13icivs.osteo[1:30,])), size = 3.5) 
volcano_cxcl13icivs.osteoplot

#save high quality plot 
ggsave("volcano_cxcl13cd4.osteovs.ici.pdf", 
       plot = volcano_cxcl13icivs.osteoplot, 
       height=4, width = 4.9, units = "in", dpi = 300)


#' figure 4B =======================================================================================================
pd1hicxcl13hicd4_icivsosteo$gene<- rownames(pd1hicxcl13hicd4_icivsosteo)
pd1hicxcl13hicd4_icivsosteo<- pd1hicxcl13hicd4_icivsosteo %>% arrange(desc(avg_log2FC))
fold_changes<- pd1hicxcl13hicd4_icivsosteo$avg_log2FC
names(fold_changes)<- pd1hicxcl13hicd4_icivsosteo$gene
Reactome <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
hallmark <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/h.all.v2023.1.Hs.symbols.gmt")
KEGG <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
GO <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c5.go.bp.v2023.1.Hs.symbols.gmt")

#combine the four gsea
Four_gene_sets <- c(hallmark, KEGG, GO, Reactome)

#getset for all pathway 
Allgenesets <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/msigdb.v2023.1.Hs.symbols.gmt")

# GSEA analysis
gsea_Myo_cxcl13icivs.oa <- fgsea(pathways = hallmark,
                              stats = fold_changes,
                              eps = 0.0,
                              minSize=15,
                              maxSize=500)

head(gsea_Myo_cxcl13icivs.oa[order(pval), ]) 

#Enrichment of CXCL13  (ICIvs.oa)
gseaicisecondvsfirst_data<- read.csv("enrichment_cxcl13.csv", stringsAsFactors = TRUE, header = TRUE)
gseaicisecondvsfirst_data<- gseaicisecondvsfirst_data %>% filter(group== "top")
gseaicisecondvsfirst.cxcl13<- ggplot(gseaicisecondvsfirst_data, aes(reorder(Pathway, NES), NES, fill= p.adjust))  +
  theme_classic() +
  geom_bar(stat = "identity")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y  = element_blank(), 
        axis.text = element_text(face = "bold", size = 10), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10), 
        legend.title = element_text(size= 12, face="bold"), 
        legend.text = element_text(size= 12, face="bold")) + 
  labs(x="Hallmark pathways", y="normalized enrichment score") + 
  coord_flip() + 
  scale_fill_gradient(low = "red", high = "blue", 
                      limits = c(0, max(gseaicisecondvsfirst_data$p.adjust, na.rm = TRUE)))
gseaicisecondvsfirst.cxcl13
ggsave("gseasecondvsfirst.cxcl13l.pdf", plot = gseaicisecondvsfirst.cxcl13, height=4, width =6, units = "in", dpi = 300)

#' figure 4C =======================================================================================================
prepsct<- PrepSCTFindMarkers(T.annotated.subsets)
CXCL13secondvsfirst<-  FindMarkers(T.annotated.subsets,assay = "SCT", ident.1 = "Second arthritis_PD-1hi CXCL13hi CD4", 
                                   ident.2 = "First arthritis_PD-1hi CXCL13hi CD4", 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25, 
                                   verbose = FALSE, 
                                   recorrect_umi = FALSE)
head(CXCL13secondvsfirst, n = 20)


# Order the Th1 data frame by adjusted p-value
volcano_cxcl13cd4 <- CXCL13secondvsfirst[order(CXCL13secondvsfirst$p_val_adj),]

# Create data frame with cut_off based on significance and avg_log2FC thresholds
results_cxcl13cd4 <- as.data.frame(mutate(as.data.frame(volcano_cxcl13cd4), 
                                    cut_off = case_when(
                                      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Downregulated",
                                      p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Upregulated",
                                      TRUE ~ "Not significant")), 
                             row.names = rownames(volcano_cxcl13cd4))

# Ensure the cut_off column is a factor with the correct levels
results_cxcl13cd4$cut_off <- factor(results_cxcl13cd4$cut_off, levels = c('Downregulated', 'Upregulated', 'Not significant'))

# Load necessary library
library(ggrepel)
genes_to_label <- c("TNF", "IFNG", "JUN", "CD69", "KLF6", "JUNB", "IER2", "NFKB1", "DNAJA1", "IRF4", "DNAJB1", "CCNL1",
                    "CXCR6", "MT1E", "AC131971.1", "TSC22D3", "MT2A")  # Replace with your gene names
volcano_cxcl13cd4 <- ggplot(results_cxcl13cd4, aes(avg_log2FC, -log10(p_val_adj))) + 
  theme_classic() +
  geom_point(aes(col = cut_off), size= 2) +
  scale_color_manual(values = c("#00AFBB","#bb0c00", "grey")) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') +
  theme(plot.title = element_blank(), 
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(face = "bold", size = 13), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 6.5), 
        legend.position = c(0.88, 0.25)) + 
  geom_text_repel(data = subset(results_cxcl13cd4, rownames(results_cxcl13cd4) %in% genes_to_label), 
                  aes(label = rownames(subset(results_cxcl13cd4,
                                              rownames(results_cxcl13cd4) %in% genes_to_label))), size = 3.5)
ggsave("volcano_cxcl13cd4.second.vsfirst.pdf", 
       plot = volcano_cxcl13cd4, 
       height=4, width = 4.9, units = "in", dpi = 300)

#' figure 4D =======================================================================================================
CXCL13secondvsfirst$gene<- rownames(CXCL13secondvsfirst)
CXCL13secondvsfirst<- CXCL13secondvsfirst %>% arrange(desc(avg_log2FC))


fold_changes<- CXCL13secondvsfirst$avg_log2FC
names(fold_changes)<- CXCL13secondvsfirst$gene
Reactome <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
hallmark <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/h.all.v2023.1.Hs.symbols.gmt")
KEGG <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
GO <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c5.go.bp.v2023.1.Hs.symbols.gmt")

#combine the four gsea
Four_gene_sets <- c(hallmark, KEGG, GO, Reactome)

#getset for all pathway 
Allgenesets <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/msigdb.v2023.1.Hs.symbols.gmt")

# GSEA analysis
gsea_Myo_cxcl13svs.f <- fgsea(pathways = hallmark,
                              stats = fold_changes,
                              eps = 0.0,
                              minSize=15,
                              maxSize=500)

head(gsea_Myo_cxcl13svs.f[order(pval), ]) 
write.xlsx(gsea_Myo_cxcl13svs.f, "gsea.cxcl13.xlsx")
gseaicisecondvsfirst_data<- read.csv("enrichment_cxcl13.csv", stringsAsFactors = TRUE, header = TRUE)
gseaicisecondvsfirst_data<- gseaicisecondvsfirst_data %>% filter(group== "top")
gseaicisecondvsfirst.cxcl13<- ggplot(gseaicisecondvsfirst_data, aes(reorder(Pathway, NES), NES, fill= p.adjust))  +
  theme_classic() +
  geom_bar(stat = "identity")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y  = element_blank(), 
        axis.text = element_text(face = "bold", size = 10), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10), 
        legend.title = element_text(size= 12, face="bold"), 
        legend.text = element_text(size= 12, face="bold")) + 
  labs(x="Hallmark pathways", y="normalized enrichment score") + 
  coord_flip() + 
  scale_fill_gradient(low = "red", high = "blue", 
                      limits = c(0, max(gseaicisecondvsfirst_data$p.adjust, na.rm = TRUE)))
gseaicisecondvsfirst.cxcl13
ggsave("gseasecondvsfirst.cxcl13l.pdf", plot = gseaicisecondvsfirst.cxcl13, height=4, width =6, units = "in", dpi = 300)

#' figure 4E =======================================================================================================
pd1hicxcl13hicd4<- ggplot(firstvssecond, aes(group,PD1hiCXCL13hiCD4_of_live_fromTreg, 
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
  labs(x = "Groups", y = "% of live cell", title = "PD1hi CXCL13hi CD4") + 
  ylim(c(0, 6))
pd1hicxcl13hicd4
ggsave("pd1hicxcl13hicd4.pdf", plot = pd1hicxcl13hicd4, height=2, width = 2, units = "in", dpi = 300)

tfh.pairedtest<- firstvssecond %>% 
  select(PD1hiCXCL13hiCD4_of_live_fromTreg, group) %>% 
  mutate(pair_id = rep(1:(nrow(.)/2), each = 2))  # Assign pair IDs   

# Reshape to wide format
df_wide <- tfh.pairedtest %>%
  pivot_wider(names_from = group, values_from = PD1hiCXCL13hiCD4_of_live_fromTreg)

# Perform paired t-test
paired_t_test <- t.test(df_wide$`first RA`, df_wide$`second RA`, paired = TRUE)

# Print the results
print(paired_t_test)  #not sure I can trust the result !! try conventional approach !! 

#' figure 4F =======================================================================================================
PD1hiCXCL13hiCD4.IFNy<- ggplot(firstvssecond, aes(group,PD1hiCXCL13hiCD4_IFNg,  
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
  labs(x = "Groups", y = "% of PD1hi CXCL13hi CD4", title = "IFN-G") + 
  ylim(c(10, 80))
PD1hiCXCL13hiCD4.IFNy
ggsave("cd8ifng.pdf", plot = PD1hiCXCL13hiCD4.IFNy, height=2, width = 2, units = "in", dpi = 300)


#' figure 4G =======================================================================================================
PD1hiCXCL13hiCD4.TNFa<- ggplot(firstvssecond, aes(group,PD1hiCXCL13hiCD4_TNFa, 
                                    shape = group, fill = group)) +
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
  theme(legend.position = "none") +
  labs(x = "Groups", y = "% of PD1hi CXCL13hi CD4", title = "TNF-A") + 
  ylim(c(10, 100))
PD1hiCXCL13hiCD4.TNFa

#' figure 4H =======================================================================================================
PD1hiCXCL13hiCD4.IL21<- ggplot(firstvssecond, aes(group,PD1hiCXCL13hiCD4_IL21, 
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
  theme(legend.position = "none") + 
  labs(x = "Groups", y = "% PD1hi CXCL13hi CD4", title = "IL-21")
PD1hiCXCL13hiCD4.IL21
ggsave("cd8il21.pdf", plot = tfh.IL21, height=2, width = 2, units = "in", dpi = 300)


#' figure 4I =======================================================================================================
sub.T.obj<- read_rds("umap.obj.pc15.rds")
DimPlot(sub.T.obj, label = TRUE) + 
  theme_minimal() + 
  NoLegend()
view(sub.T.obj@meta.data)
FeaturePlot(sub.T.obj, features = c("CD4", "CD8A"))

# Subset CD4 
seurat.obj.CD4<- subset(sub.T.obj, idents = c("0", "4", "6", "8", "10"))
DimPlot(seurat.obj.CD4, reduction= "umap", label = TRUE) +
  theme_classic() +
  NoLegend()

#remove contaminated cells! 
plot<- DimPlot(seurat.obj.CD4, reduction= "umap", label = TRUE) +
  theme_minimal() +
  NoLegend()
seurat.object<- CellSelector(plot = plot, seurat.obj.CD4)

DimPlot(seurat.object, reduction= "umap", label = TRUE) +
  theme_minimal() +
  NoLegend()

#subset selectedcell !
seurat.obj.CD4 <- subset(seurat.object, idents = c("0", "4", "6", "8", "10", "15"))
DimPlot(seurat.obj.CD4, reduction= "umap", label = TRUE) +
  theme_minimal() +
  NoLegend()

# SCT/PCA 
cd4.subcluster<- SCTransform(seurat.obj.CD4, 
                         method= "glmGamPoi")
cd4.subcluster<- RunPCA(cd4.subcluster, verbose = FALSE)
saveRDS(cd4.subcluster, file = "cd4.subcluster.sct.pca.rds")  #save on 01/31/2024
T.subcluster<- read_rds("T.subcluster.sct.pca.rds")

ElbowPlot(cd4.subcluster, ndims = 50,reduction = "pca") +
  ggtitle("Number of principle component")  + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Integration/clustering
cd4.T.integration<- RunHarmony(cd4.subcluster,
                              group.by.vars = "orig.ident")
saveRDS(cd4.T.integration, file = "cd4.cell.harmony.rds")   #save on 01/31/2024

CD4.T.pc15 <- FindNeighbors(cd4.T.integration, 
                                 reduction = "harmony", 
                               dims = 1:15)

CD4.T.pc15<- FindClusters(CD4.T.pc15 ,
                                 resolution = 0.205)
CD4.T.pc15<- RunUMAP(object = CD4.T.pc15,
                          reduction = "harmony",
                          dims = 1:15)
CD4_subclu<- DimPlot(CD4.T.pc15, 
                     reduction= "umap", 
                     label = FALSE, 
                     sizes.highlight = 10) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        panel.grid = element_blank(),
        axis.text = element_text(face = "bold", size= 16), 
        axis.title = element_text(face = "bold", size= 18), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) + 
  NoLegend()
CD4_subclu
ggsave("CD4.pdf", plot = CD4_subclu, height=3.5, width = 4, units = "in", dpi = 300)

#' figure 4I =======================================================================================================
cd4.monocle<- readRDS("CD4.obj.rds")
cd4_subclu<- DimPlot(cd4.monocle, 
                    reduction= "umap", 
                    label = FALSE, 
                    sizes.highlight = 20, 
                    pt.size = 1.5) +
  theme_minimal() +
  NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(face = "bold", size= 16), 
        axis.title = element_text(face = "bold", size= 18),
         panel.grid.major = element_line(size = 1), 
        panel.grid.minor = element_line(size= 1))
cd4_subclu

#Second arthritis 
second_arthritis.cd4<- subset(cd4.monocle, subset = type == "Second arthritis")
monocle.cds.cd4s <- as.cell_data_set(second_arthritis.cd4)
monocle.cds.cd4s <- cluster_cells(cds = monocle.cds.cd4s, 
                             reduction_method = "UMAP")
monocle.cds.cd4s <- learn_graph(monocle.cds.cd4s)
startCells <- Cells(second_arthritis.cd4)[second_arthritis.cd4$seurat_clusters == 3]
monocle.cds.cd4s <- order_cells(monocle.cds.cd4s,
                           reduction_method = "UMAP", 
                           root_cells = startCells)
plot_cells(cds = monocle.cds.cd4s,
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

monocle.cds.cd4s$monocle3_pseudotime <- pseudotime(monocle.cds.cd4s)
data.pseudo_second <- as.data.frame(colData(monocle.cds.cd4s))
ggplot(data.pseudo_second , aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()
pseudotime_second<- ggplot(data.pseudo_second, aes(x = monocle3_pseudotime, y = reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) +
  geom_density_ridges() + 
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title.y  = element_blank(), 
        axis.text = element_text(size= 15, face = "bold"), 
        axis.title.x = element_text(siz= 18, face = "bold"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.3), 
        panel.grid = element_blank()) + 
  labs(x= "Pseudotime")
pseudotime_second
ggsave("cd4.pseudo.s.pdf", plot = pseudotime_second, height=4.5, width = 4, units = "in", dpi = 300)

#First arthritis pseudotime 
first_arthritis.cd4<- subset(cd4.monocle, subset = type == "First arthritis")
monocle.cds.cd4f <- as.cell_data_set(first_arthritis.cd4)
monocle.cds.cd4f <- cluster_cells(cds = monocle.cds.cd4f, 
                             reduction_method = "UMAP")
monocle.cds.cd4f <- learn_graph(monocle.cds.cd4f)
startCells <- Cells(first_arthritis.cd4)[first_arthritis.cd4$seurat_clusters == 3]
monocle.cds.cd4f<- order_cells(monocle.cds.cd4f,
                           reduction_method = "UMAP", 
                           root_cells = startCells)
plot_cells(cds = monocle.cds.cd4f,
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

monocle.cds.cd4f$monocle3_pseudotime<- pseudotime(monocle.cds.cd4f)
data.pseudo_first <- as.data.frame(colData(monocle.cds.cd4f))
pseudotime_first<- ggplot(data.pseudo_first, aes(x = monocle3_pseudotime, y = reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) +
  geom_density_ridges() + 
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title.y  = element_blank(), 
        axis.text = element_text(size= 15, face = "bold"), 
        axis.title.x = element_text(siz= 18, face = "bold"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.3), 
        panel.grid = element_blank()) + 
  labs(x= "Pseudotime")
pseudotime_first
ggsave("cd4.pseudo.f.pdf", plot = pseudotime_first, height=4.5, width = 4, units = "in", dpi = 300)

#' End-of-figure 4 ===================================================================================================================
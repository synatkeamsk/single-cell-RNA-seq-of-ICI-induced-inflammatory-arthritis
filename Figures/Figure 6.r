
require(tidyverse)
library(Seurat)
library(ggrepel)
library(fgsea)
library(patchwork)
Treg.obj<- readRDS("Treg.obj.annotated_one.rds")

#subset the first and second IA !!! 
Treg.obj<- subset(Treg.obj, subset = type %in% c("First arthritis","Second arthritis")) 

#' Figure 6A   ==================================================================================================
Tregone<- DimPlot(Treg.obj, label = FALSE, 
                  group.by = "type", 
                  pt.size = 0.1) + 
  theme_classic()  +
    theme(plot.title = element_blank(),
        axis.text = element_text(face = "bold", size= 11), 
        axis.title = element_blank(), 
        legend.position = "none") 
Tregone

#' Figure 6B  ====================================================================================================
Treg_frequency<- read.csv("cluster_frequency/cluster_frequency.csv", stringsAsFactors = TRUE)
Treg_frequency<- Treg_frequency %>% 
  filter(Arthritis %in% c("first arthritis", "second arthritis")) %>% 
  filter(cluster_name == "Treg")

Tregfreq<- ggplot(Treg_frequency, aes(x= Arthritis, y=frequency, fill= Arthritis)) + 
  theme_classic() +
  geom_violin(trim = FALSE) +
  geom_point(size= 5, aes(x=Arthritis), shape=21, position = position_dodge(width = 0)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14
  ), axis.title.y = element_text(size= 20, face = "bold"), 
  axis.title.x = element_blank(),
  axis.text = element_text(size=16, face = "bold"), 
  axis.text.x = element_blank()) + 
  theme(legend.position = "none") + 
  ylab("% of total cells") + 
  ggtitle("Regulatory T cell")
Tregfreq

#' Figure 6C ====================================================================================================
features<- c( "CTLA4", "PDCD1", "TIGIT","LAG3","IL2RA", "ICOS","KLRG1", 
             "IFNG", "TNF", "IL1B", "IL17", "IL10", "TGFB1", "IKZF2", "FOXP3","CD4")


dotplot<- DotPlot(Treg.T.pc10_filter, 
                         features = features, 
                         cols = c("blue", "red"),
                         dot.scale = 6, 
                         group.by = "type") + 
  theme_classic() +
  RotatedAxis() + 
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15), 
        axis.title.y= element_blank(), 
        axis.text = element_text(size=13, face = "bold"), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.title = element_text(size= 14, face = "bold"), 
        legend.text = element_text(size= 14, face = "bold")) + 
  theme(legend.position = "right")
ggsave("allTreg.pdf", plot = dotplot, height=3, width = 6,  units = "in", dpi = 300)

#'  Figure 6D   ===================================================================================================
Treg.obj$celltypes<-Idents(Treg.obj)
Treg.obj$arthritis_celltype <- paste(Treg.obj$type, sep = "_", Idents(Treg.obj))
Idents(Treg.obj) <- "arthritis_celltype"
prepsct<- PrepSCTFindMarkers(Treg.obj)
Tregsecondvsfirst<-  FindMarkers(prepsct,assay = "SCT", ident.1 = "Second arthritis_Treg", 
                                   ident.2 = "First arthritis_Treg", 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25, 
                                   verbose = FALSE, 
                                   recorrect_umi = FALSE)
Tregorder<- Tregsecondvsfirst %>% arrange(desc(avg_log2FC))
Tregorder<-Tregorder[Tregorder$p_val_adj<0.05, ]
write.csv(Tregorder, file = "deg_tregall.csv")
# Order the Th1 data frame by adjusted p-value
volcano_Treg <- Tregsecondvsfirst[order(Tregsecondvsfirst$p_val_adj),]

# Create data frame with cut_off based on significance and avg_log2FC thresholds
results_Treg <- as.data.frame(mutate(as.data.frame(volcano_Treg), 
                                    cut_off = case_when(
                                      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Downregulated",
                                      p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Upregulated",
                                      TRUE ~ "Not significant")), 
                             row.names = rownames(volcano_Treg))

# Ensure the cut_off column is a factor with the correct levels
results_Treg$cut_off <- factor(results_Treg$cut_off, levels = c('Downregulated', 'Upregulated', 'Not significant'))

# Load necessary library
options(ggrepel.max.overlaps = Inf)
genes_to_label <- c("TNF", "IFNG", "JUN", "CD69", "KLF6", "JUNB", "IER2", "NFKB1", "DNAJA1", "IRF4", "DNAJB1", "CCNL1",
                    "CXCR6", "MT1E", "AC131971.1", "TSC22D3", "MT2A")  # Replace with your gene names

volcano_Treg_plot <- ggplot(results_Treg, aes(avg_log2FC, -log10(p_val_adj))) + 
  theme_classic() +
  geom_point(aes(col = cut_off), size= 2) +
  scale_color_manual(values = c("#00AFBB","#bb0c00", "grey")) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') +
  theme(plot.title = element_blank(), 
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(face = "bold", size = 13), 
        legend.title = element_blank(), 
        legend.text = element_blank()) + 
  geom_text_repel(data = results_Treg[1:25,], 
                  aes(label = rownames(results_Treg[1:25,])), size = 2) + 
   theme(legend.position = "none") +
  ggtitle("Second vs. first arthritis (PD1hi CXCL13hi CD4)")
volcano_Treg_plot
ggsave("volcano.Tregall.pdf", plot = volcano_Treg_plot, height=3, width = 3.5, units = "in", dpi = 300)

#'  Figure 6E   ===================================================================================================
Tregsecondvsfirst$gene<- rownames(Tregsecondvsfirst)
Tregsecondvsfirst<- Tregsecondvsfirst %>% arrange(desc(avg_log2FC))
fold_changes<- Tregsecondvsfirst$avg_log2FC
names(fold_changes)<- Tregsecondvsfirst$gene
Reactome <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
hallmark <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/h.all.v2023.1.Hs.symbols.gmt")
KEGG <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
GO <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/c5.go.bp.v2023.1.Hs.symbols.gmt")
Four_gene_sets <- c(hallmark, KEGG, GO, Reactome)
Allgenesets <- fgsea::gmtPathways("Y:/Recurrent__Arthritis__Single cell/Single_Cell_Arthritis/msigdb.v2023.1.Hs.symbols.gmt")

# GSEA analysis
gsea_treg <- fgsea(pathways = hallmark,
                              stats = fold_changes,
                              eps = 0.0,
                              minSize=15,
                              maxSize=500)

head(gsea_Myo_treg[order(pval), ]) 

# make an enrichment plot for a pathway:
TNF<- plotEnrichment(hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
                     fold_changes) + 
  theme_classic()

IL2<- plotEnrichment(hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]],
                     fold_changes) +  
  theme_classic()

IFNG<- plotEnrichment(hallmark[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
                     fold_changes) + 
  theme_classic()
Treg<- TNF + IL2 + IFNG
ggsave("gase_Treg.tiff", plot = Treg, height=2.5, width =5, units = "in", dpi = 300)

#'  Figure 6F  =========================================================================================================
Tregfirstvssecond<- ggplot(firstvssecond, aes(x = group, y = Treg_1_of_live)) +
  geom_bar(aes(fill = group), stat = "summary", fun = "mean", width = 0.6, alpha = 0.7) +
  geom_point(aes(shape = group, fill = group), 
             position = position_jitter(width = 0.15), 
             size = 3) +
  geom_line(aes(group = paired), color = "black") +
  scale_shape_manual(values = c(16, 17)) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "none"
  ) +
  labs(y = "% of CD4", title = "Th1")
Tregfirstvssecond
ggsave("Treg.first_second.pdf", plot = Tregfirstvssecond, height=1.5, width = 1.5, units = "in", dpi = 300)

#'  Figure 6G  =========================================================================================================
Treg.TNF<- ggplot(firstvssecond, aes(x = group, y = Treg1_TNF)) +
  geom_bar(aes(fill = group), stat = "summary", fun = "mean", width = 0.6, alpha = 0.7) +
  geom_point(aes(shape = group, fill = group), 
             position = position_jitter(width = 0.15), 
             size = 3) + 
  geom_line(aes(group = paired), color = "black") +
  
  scale_shape_manual(values = c(16, 17)) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "none"
  )
Treg.TNF
ggsave("Treg.TNF.pdf", plot = Treg.TNF, height=1.5, width = 1.5, units = "in", dpi = 300)

#'  Figure 6H  =========================================================================================================
Treg.ifng<- ggplot(firstvssecond, aes(x = group, y = Treg1_IFNG)) +

  geom_bar(aes(fill = group), stat = "summary", fun = "mean", width = 0.6, alpha = 0.7) +
  geom_point(aes(shape = group, fill = group), 
             position = position_jitter(width = 0.15), 
             size = 3) + 
  geom_line(aes(group = paired), color = "black") +
  
  scale_shape_manual(values = c(16, 17)) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_blank(),
    legend.position = "none"
  )
Treg.ifng
ggsave("Treg.ifng.pdf", plot = Treg.ifng, height=1.5, width = 1.5, units = "in", dpi = 300)

#' End-of-figure 6 ===================================================================================================================
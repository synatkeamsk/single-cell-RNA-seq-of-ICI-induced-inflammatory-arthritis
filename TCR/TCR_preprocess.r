suppressMessages(require(tidyverse))
suppressMessages(library(scRepertoire))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

#Read the data 
A164S <- read.csv("164S/filtered_contig_annotations_164S.csv")
A164S2 <- read.csv("164S2/filtered_contig_annotations_164S2.csv")
A184S2 <- read.csv("184S2/filtered_contig_annotations_184S2.csv")
A184S3 <- read.csv("184S3/filtered_contig_annotations_184S3.csv")
A218S <- read.csv("218S/filtered_contig_annotations_218S.csv")
A218S2 <- read.csv("218S2/filtered_contig_annotations_218S2.csv")
ASA5 <- read.csv("SA5/filtered_contig_annotations_SA5.csv")
ASA6 <- read.csv("SA6/filtered_contig_annotations_SA6.csv")
ASA7 <- read.csv("SA7/filtered_contig_annotations_SA7.csv")

#Contig_list
contig_list <- list(A164S, A164S2, A184S2, A184S3, A218S, A218S2, ASA5, ASA6, ASA7)

#Combine TCR
combined <- combineTCR(contig_list, 
                       samples = c("164S", "164S2", "184S2", "184S3", "218S","218S2", "SA5", "SA6", "SA7"),
                       removeNA = TRUE, 
                       removeMulti = TRUE,
                       filterMulti = TRUE)

combined <- addVariable(combined, variable.name =  "patientid",
                       variables = c("p164S", "p164S2", "p184S2", "p184S3", "p218S","p218S2", "pSA5", "pSA6", "pSA7"))

combined <- addVariable(combined, variable.name = "type",
                       variables = c("First arthritis", "Second arthritis", "First arthritis", "Second arthritis", 
                                     "First arthritis", "Second arthritis", "OA", "OA", "OA"))

combined <- addVariable(combined, variable.name = "group",
                       variables = c("ICI-arthritis", "ICI-arthritis", "ICI-arthritis", "ICI-arthritis", 
                                     "ICI-arthritis", "ICI-arthritis", "OA", "OA", "OA"))

                                    
#Integration with ScRNA seq ! 
library(Seurat)
seurat.obj.T<- readRDS("umap.obj.pc15.rds")

#Dimplot 
DimPlot(seurat.obj.T)

#create new column variable! 
seurat.obj.T@meta.data$group <- ifelse(seurat.obj.T@meta.data$type %in% c("First arthritis", "Second arthritis"), "ICI-arthritis", "Osteoarthritis")
view(seurat.obj.T@meta.data)

#modify the clusters
seurat.obj.T<- RenameIdents(seurat.obj.T, `0` ="c0",`10`= "c1", `4`= "c2", 
                                   `6`= "c3", `8`= "c4", `15`= "c5", `18`= "c5", `3` ="c6", `11`="c6", 
                                  `2` ="c7", `7`= "c7", `5`= "c8", `1` ="c9", `9`= "c9",
                                   `14`= "c10", `17`= "c11", `13`= "c12", `12`= "c13", `16`= "c14") 

cell.barcodes.T <- rownames(seurat.obj.T[[]])
cell.barcodes.T <- stringr::str_split(cell.barcodes.T, "_", simplify = TRUE)[,1] 
cell.barcodes.T <- paste0(seurat.obj.T$orig.ident, "_", cell.barcodes.T) 
seurat.obj.T <- RenameCells(seurat.obj.T, new.names = cell.barcodes.T)

#Combine Expression 
integrated.obj.T <- combineExpression(combined, seurat.obj.T, 
                            cloneCall = "gene", 
                            proportion = FALSE, 
                            cloneSize  = c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

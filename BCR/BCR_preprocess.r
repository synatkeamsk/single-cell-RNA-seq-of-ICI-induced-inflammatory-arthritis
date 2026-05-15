suppressMessages(require(tidyverse))
suppressMessages(library(scRepertoire))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

#Read the data
A184S2.bcr <- read.csv("184S2/filtered_contig_annotations_184S2_bcr.csv")
A184S3.bcr <- read.csv("184S3/filtered_contig_annotations_184S3_bcr.csv")
A218S.bcr <- read.csv("218S/filtered_contig_annotations_218S_bcr.csv")
A218S2.bcr <- read.csv("218S2/filtered_contig_annotations_218S2_bcr.csv")

#create contig_list
contig_list.bcr <- list(A184S2.bcr, A184S3.bcr, A218S.bcr, A218S2.bcr)

#Combine BCR
combined.bcr<- combineBCR(contig_list.bcr, 
                       samples = c("184S2", "184S3", "218S","218S2"),
                       removeNA = TRUE, 
                       removeMulti = TRUE,
                       filterMulti = TRUE)

#Add variables 
combined.bcr <- addVariable(combined.bcr, variable.name =  "patientid",
                       variables = c("p184S2", "p184S3", "p218S","p218S2"))

combined.bcr <- addVariable(combined.bcr, variable.name = "type",
                       variables = c("First arthritis", "Second arthritis", 
                                     "First arthritis", "Second arthritis"))

combined.bcr <- addVariable(combined.bcr, variable.name = "group",
                       variables = c("ICI-arthritis", "ICI-arthritis", "ICI-arthritis", "ICI-arthritis"))       

# Integration BCR with scRNAseq
seurat.obj.B<- readRDS("B.cell.rds")
DimPlot(seurat.obj.B, label = TRUE)
seurat.obj.B<- subset(seurat.obj.B, subset = type %in% c("First arthritis","Second arthritis")) 
cell.barcodes.B <- rownames(seurat.obj.B[[]])
cell.barcodes.B <- stringr::str_split(cell.barcodes.B, "_", simplify = TRUE)[,1] 
cell.barcodes.B <- paste0(seurat.obj.B$orig.ident, "_", cell.barcodes.B) 
seurat.obj.B <- RenameCells(seurat.obj.B, new.names = cell.barcodes.B)


# Integration of B cell cluster and BCR
integrated.obj.B <- combineExpression(combined.bcr, seurat.obj.B, 
                                        cloneCall = "gene", 
                                        proportion = FALSE, 
                                        cloneSize  = c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

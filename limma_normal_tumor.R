source("helper_functions.R")

limma_paired <- function(my_IDs,dataframe,limma_name,up_name,down_name){
  
  condition <- as.factor(my_IDs$condition)
  patientID <- as.factor(my_IDs$participant)
  
  #design matrix
  design.matrix <- model.matrix(~0+condition+patientID)
  colnames(design.matrix)[c(1,2)] <- c("cancer","normal")
  
  #voom transformation
  dataframe <- voom(dataframe,design.matrix,plot=TRUE)
  
  # Making group contrasts 
  N_C_cotr <- makeContrasts("cancer-normal", levels= design.matrix)
  
  # Filter for significance - set log fold change and fdr cutoff
  N_C_L <- DE_limma(N_C_cotr, dataframe, design.matrix, 1, 0.01)
  
  # differentially expressed genes file
  write.csv(N_C_L, limma_name, quote = FALSE)
  
  #number of up-regulated genes cancer vs normal
  print(paste0("the analysis detected ",length(which(N_C_L$direction == "up"))," up-regulated genes"))
  #number of down-regulated genes cancer vs normal
  print(paste0("the analysis detected ",length(which(N_C_L$direction== "down")), " down-regulated genes"))
  
  # up and down regulated genes
  up <- data.frame(rownames(N_C_L[N_C_L$direction == "up", ]))
  down <- data.frame(rownames(N_C_L[N_C_L$direction == "down", ]))
  
  write.table(up, up_name, sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  write.table(down, down_name, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
}

#--------------------------------------------------------------------------
# LUAD paired after removing low tumor purity samples
#-------------------------------------------------------------------------

dataframe_LUAD <- get(load("LUAD/paired/LUAD_PreprocessedData_paired_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
#check the tumor samples number
length(which(my_IDs$condition=="cancer"))
#check the normal samples number
length(which(my_IDs$condition=="normal"))
limma_name <- "limma_LUAD_paired_tumorPurity.csv"
up_name <- "limma/up_limma_LUAD_paired_tumorPurity.txt"
down_name <- "limma/down_limma_LUAD_paired_tumorPurity.txt"

limma_paired(my_IDs,dataframe_LUAD,limma_name,up_name,down_name)

#----------------------------------------------------------------------------------
# Create Scatterplots
#----------------------------------------------------------------------------------
up_LUAD <- read.csv("limma_LUAD_paired_tumorPurity.csv", row.names = 1)
up_LUAD$symbol <- rownames(up_LUAD)
#scatterplot to compare logFC
TCGAVisualize_volcano(up_LUAD$logFC, up_LUAD$adj.P.Val,
                      filename = "LUAD_limma.pdf", xlab = "logFC",
                      names = up_LUAD$symbol, show.names = "highlighted",
                      x.cut = 1, y.cut = 0.05, 
                      highlight = (up_LUAD$symbol)[which(abs(up_LUAD$logFC) >= 2)],
                      highlight.color = "red"
                      )

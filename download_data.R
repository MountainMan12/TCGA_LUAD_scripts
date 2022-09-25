#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks", "SummarizedExperiment", "EDASeq")


library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)

#---------------------------------------------------
# download all samples - LUAD
#---------------------------------------------------

query_LUAD <- GDCquery(project = "TCGA-LUAD",
                       data.category = "Gene expression",
                       data.type = "Gene expression quantification",
                       platform = "Illumina HiSeq",
                       file.type = "results",
                       legacy = TRUE,
                       experimental.strategy = "RNA-Seq",
                       sample.type = c("Primary Tumor","Solid Tissue Normal"))

GDCdownload(query_LUAD)

data_LUAD <- GDCprepare(query_LUAD, save = TRUE, save.filename = "./LUAD/all/LUAD_Illumina_HiSeq_all.rda")

# Which samples are primary solid tumor
dataSmTP <- TCGAquery_SampleTypes(query_LUAD$results[[1]]$cases,"TP") 
# which samples are solid tissue normal
dataSmNT <- TCGAquery_SampleTypes(query_LUAD$results[[1]]$cases,"NT")

matched <- TCGAquery_MatchedCoupledSampleTypes(c(dataSmTP, dataSmNT), c("TP","NT"))

query_LUADpaired <- GDCquery(project = "TCGA-LUAD", 
                              legacy = TRUE,
                              data.category = "Gene expression",
                              data.type = "Gene expression quantification",
                              platform = "Illumina HiSeq", 
                              file.type = "results",
                              experimental.strategy = "RNA-Seq",
                              sample.type = c("Primary Tumor","Solid Tissue Normal"),
                              barcode = matched)

GDCdownload(query_LUADpaired)

data_LUADpaired <- GDCprepare(query = query_LUADpaired, save = TRUE, save.filename = "./LUAD/paired/LUAD_Illumina_HiSeq_paired.rda")

length(which(colData(data_LUADpaired)$shortLetterCode =="TP"))
length(which(colData(data_LUADpaired)$shortLetterCode =="NT"))

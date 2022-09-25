#Get all clinical data for paired samples
clinDA <- as.data.frame(data_LUADpaired@colData)
#Filter only samples in the gene expression data 
clinDa <- clinDA[colnames(dataframe_LUAD),]

TCGAanalyze_survival(
  data = clinDa,
  clusterCol = "ajcc_pathologic_stage",
  main = "TCGA Set\n LUAD",
  conf.int = F,
  height = 10,
  width=10
)


#COX-PH Survival Analysis
library(stringr)

rnaExpr_gene <- dataframe_LUAD[up_LUAD$symbol[1:10],]
rnaExpr_gene_norm <- log2(rnaExpr_gene)+1

tabSurvKM <- TCGAanalyze_SurvivalKM(clinical_patient = clinDa,
                                    dataGE = rnaExpr_gene_norm,
                                    Genelist = rownames(rnaExpr_gene_norm),
                                    Survresult = TRUE,
                                    ThreshTop = 0.67,
                                    p.cut = 0.9,
                                    ThreshDown = 0.33)
  


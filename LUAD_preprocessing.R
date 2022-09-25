SE_LUAD<-get(load("./LUAD/paired/LUAD_Illumina_HiSeq_paired.rda"))

#check the samples number
length(which(colData(SE_LUAD)$shortLetterCode =="TP"))
length(which(colData(SE_LUAD)$shortLetterCode =="NT"))

#remove samples with low tumor purity
list <- TCGAtumor_purity(colnames(SE_LUAD),0,0,0,0,0.6)

# Get number of pure barcodes
length(list$pure_barcodes)

# Get number of barcodes to be filtered out
length(list$filtered)

SE_LUAD <- SE_LUAD[,colnames(SE_LUAD) %in% union(list$pure_barcodes,list$filtered)]
#get paired samples only
paired <- TCGAquery_MatchedCoupledSampleTypes(colnames(SE_LUAD),c("NT","TP"))
SE_LUAD <- SE_LUAD[,colnames(SE_LUAD) %in% paired]
save(SE_LUAD, file = "./LUAD/paired/LUAD_Illumina_HiSeq_paired_tumorPurity.rda")

#Preprocessing
dataPrep <- TCGAanalyze_Preprocessing(object = SE_LUAD, cor.cut = 0.6)

#-------Quality control--------------
#######The dataset is usually filtered at this stage to remove any######### 
#######genes that are not expressed. It can reduce######################### 
#######the time and memory required to perform some of the analysis########
#######Let’s say that for a gene to be “expressed” in a particular sample## 
#######################we need to see 5 or more counts#####################

is_expressed <- dataPrep >= 5

df <- data.frame(Expressed = rowSums(is_expressed))
ggplot(df, aes(x=Expressed)) + geom_bar()

# Keep genes that are expressed in atleast 5 samples
keep <- rowSums(dataPrep >= 5) >= 5
table(keep)

#Remove genes that do not fit the inclusion criteria
dataPrep <- dataPrep[keep,]

#Visualise the count distributions
boxplot(dataPrep)

#For easy interpretation we take log10 of all values
boxplot(log10(dataPrep))

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

save(dataFilt, file = "./LUAD/paired/LUAD_PreprocessedData_paired_tumorPurity.rda")

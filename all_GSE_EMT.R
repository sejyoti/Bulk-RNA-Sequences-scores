
library(readxl)
library(matlabr)
#library(xlsx)
library(data.table)
library(ggplot2)
getDTthreads()
setDTthreads(4)

source("/home/sejyoti/EMT_Scoring_RNASeq/UsefulFunctions.R")
# source("/media/csb/New Volume/joel/FPKMtoTPM.v2.R")
source("/home/sejyoti/EMT_Scoring_RNASeq/EMT_Scoring_RNASeq-master/counts_to_TPM.R")
source("/home/sejyoti/EMT_Scoring_RNASeq/EMT_Scoring_RNASeq-master/EMT_score_func.R")

## Folder with all the raw read counts
setwd("/home/sejyoti/EMT_Scoring_RNASeq/EMT_Scoring_RNASeq-master")
GSE = "GSE85857"

setwd(paste("/home/sejyoti/EMT_Scoring_RNASeq/EMT_Scoring_RNASeq-master/Data/", GSE, "/", sep = ""))

fileList = list.files(pattern = ".txt")

counts = lapply(fileList, read.delim, header = T, sep = "\t") #change delimiter "comma" or "tab" accordingly
gseIDs = sapply(strsplit(fileList, split = "_"), function(x) x[1])
setwd("/home/sejyoti/EMT_Scoring_RNASeq/EMT_Scoring_RNASeq-master/Data/")
counts[[1]] = KeepMaxGene(counts[[1]]); counts[[1]] <- cbind(rownames(counts[[1]]), counts[[1]])

for(dataNum in 1:length(counts)){
  #-------------------------------------------------------
  #if TPM given
  log2_TPM <- log2(data.frame(counts[[1]], row.names = 1)+1)
  fwrite(log2_TPM, paste("../Data_generated/", GSE, "_TPM.tsv", sep = ""), sep = "\t", row.names = T)
  #------------------------
  MA_val = rnaToMA(log2_TPM)
  EMT76GS_Score = EMT76GS(MA_val)
  KS_score=KSScore(MA_val)
  MLR_score = mlrEMTPred(MA_val,gseIDs[dataNum])
}

# 76GS and KS_score <- fread(paste("../Data_generated/", GSE, "_emt_score_GS76.txt", sep = ""), header = T, sep = "\t")

EMT_scores_all <- data.frame("Sample" = EMT76GS_Score[,1], 
                             "GS76" = as.numeric(EMT76GS_Score[,2]), 
                             "KS" = KS_score[,1])


fwrite(EMT_scores_all, paste("../Output/", GSE, "_all_score.txt", sep = ""), sep = "\t")






########################################################################## 
#' The main program for DECODE algorithm
#'
#' To run an example using expression data with 1400 genes.
#' 
#' runDecode("\\extdata\\geneSet.txt","\\extdata\\Expression_data_1400genes.txt")
#'
#' or 
#'
#' runDecode("/extdata/geneSet.txt","/extdata/Expression_data_1400genes.txt")
#'
#' The sample data with 1400 genes takes 16 minutes to complete. (Computer used: An Intel Core i7-4600 processor, 2.69 GHz, 8 GB RAM)
#' @title Differential Co-Expression and Differential Expression Analysis
#' @description Given a set of gene expression data and functional gene set data, the program will return a table summary for the selected gene sets with high differential co-expression and high differential expression (HDC-HDE). User need to specify the input paths for the gene expression data and functional gene set data.
#' @param geneSetInputFile Path for functional gene set data
#' @param geneExpressionFile Path for gene expression data
#'
#'
#'
#' Input: 
#' (1) gene expression data
#'
#' (2) functional gene set data
#'
#' Output: Table summary for the selected HDC-HDE gene sets, 'out_summary.txt' 
#'
#' 
#' 
#' Data format for gene expression data (Columns are tab-separated):
#'
#' 
#'  Column 1: Official gene symbol
#'
#'  Column 2: Probe ID
#'
#'  Starting from column 3: Expression for different samples
#'
#'
#'  Row 1 (starting from column 3): Sample class ("1" indicates control group; "2" indicates case group)
#'
#'  Row 2: Sample id
#'
#'  Starting from row 3: Expression for different genes  
#'
#'
#'  Example:
#'
#' geneName   probeID         2       2        2        1         1         1
#'
#' -            -          Case1    Case2    Case3   Control1  Control2  Control3
#'
#' 7A5      ILMN_1762337  5.12621  5.19419  5.06645  5.40649   5.51259   5.387
#'
#' A1BG     ILMN_2055271  5.63504  5.68533  5.66251  5.37466   5.43955   5.50973
#'
#' A1CF     ILMN_2383229  5.41543  5.58543  5.43239  5.49634   5.62685   5.36962
#'
#' A26C3    ILMN_1653355  5.56713  5.5547   5.59547  5.46895   5.49622   5.50094
#'
#' A2BP1    ILMN_1814316  5.23016  5.33808  5.31413  5.30586   5.40108   5.31855
#'
#' A2M      ILMN_1745607  7.65332  6.56431  8.20163  9.19837   9.04295   10.1448
#'
#' A2ML1    ILMN_2136495  5.53532  5.93801  5.33728  5.36676   5.79942   5.13974
#'
#' A3GALT2  ILMN_1668111  5.18578  5.35207  5.30554  5.26107   5.26536   5.28932
#'  
#' 
#'  
#'  
#' Data format for functional gene set data (Columns are tab-separated):
#'
#' 
#'  Column 1: Functional gene set name
#'
#'  Column 2: Other description such as gene set id
#'
#'  Starting from column 3: Official gene symbols for the functional gene set
#'
#'  Example:
#'
#' B cell activation   GO\\GO:0042113   AKAP17A   ZAP70   PFDN1 ...
#'
#' apoptotic signaling pathway   GO\\GO:0097190   ITPR1   PTH   DNAJC10   HINT1 ...
#'
#'
#' @import utils stats
#'@examples
#'\dontrun{
#' path = system.file('extdata', package='decode')
#' geneSetInputFile = file.path(path, "geneSet.txt")
#' geneExpressionFile = file.path(path, "Expression_data_50genes.txt")
#' runDecode(geneSetInputFile, geneExpressionFile)
#'}
#' @export
##########################################################################
# DECODE  (c) Copyright 2014 by The Hong Kong Polytechnic University, Department of Health Technology and Informatics 
# Written by Thomas Lui
# Permission is granted to copy and use this program provided no fee is
# charged for it and provided that this copyright notice is not removed.
##########################################################################
runDecode =function(geneSetInputFile, geneExpressionFile) {


  # gene set data
#  geneSetInputFile = 'inst\\extdata\\geneSet.txt'
  # gene expression data
#  geneExpressionFile ="inst\\extdata\\Expression_data.txt"
  # significant level
  significanceLevel =0.05
  
  
  ###########################
  # read expression data
  ###########################
  print("Reading gene expression data...")
  rawData = read.table(geneExpressionFile, sep="\t",row.names = NULL)
  
  # number of row
  # dimension of the data  
  dim_rawData=dim(rawData)
  
  MaxRow =nrow(rawData)
  # number of column
  MaxCol =ncol(rawData)

  # dimension of the data  
  dim_rawData=dim(rawData)
  
  # number of row
  MaxRow =nrow(rawData)
  # number of column
  MaxCol =ncol(rawData)
  # remove last column if empty
  if (is.na(rawData[1,MaxCol])) {
    rawData =rawData[,-MaxCol]
  MaxCol=MaxCol-1
  }
  
  # get gene name
  geneName=as.matrix(rawData[3:MaxRow,1])
  # get expression matrix
  expressionMatrix_both = as.matrix(rawData[3:MaxRow,3:MaxCol])
  expressionMatrix_both=matrix(as.double(expressionMatrix_both),(MaxRow-2),(MaxCol-2))
  # initialize expression matrix for normal group
  expressionMatrix_normal=matrix(nrow=(MaxRow-2),ncol=0)
  # initialize expression matrix for disease group
  expressionMatrix_disease=matrix(nrow=(MaxRow-2),ncol=0)
  
  # get expression matrix for both groups
  for (j in 3:MaxCol) {
    if (rawData[1,j]==1) {
      expressionMatrix_normal =cbind(expressionMatrix_normal,expressionMatrix_both[,(j-2)])
  }
    if (rawData[1,j]==2) {
      expressionMatrix_disease =cbind(expressionMatrix_disease,expressionMatrix_both[,(j-2)])
  }
  }
  # dimension of expression matrix
  dim_expressionMatrix_both=dim(expressionMatrix_both)
  dim_expressionMatrix_both
  # number of genes  
  MaxGene =nrow(expressionMatrix_both)
  # number of samples 
  MaxSample =ncol(expressionMatrix_both)


  ########################################### 
  # Calculate t-test and fold change
  ###########################################
  print("Calculating t-statistics...")
  foldChange = rowMeans(2^expressionMatrix_disease,na.rm=TRUE) / rowMeans(2^expressionMatrix_normal,na.rm=TRUE)
  t_result=data.frame(tScore=double(),pValue=double(),absTScore=double())
  for (i in 1:MaxGene) {
    if (i %% 1000==0) {
      #print(i)
    }
    # ttest - assume unequal variance
    temp_t_result=t.test(expressionMatrix_disease[i,],expressionMatrix_normal[i,])  #automaticcally ignore missing values 
    t_result[i,"tScore"] = as.double(temp_t_result$statistic)
    t_result[i,"pValue"] = as.double(temp_t_result$p.value)
    t_result[i,"absTScore"] = abs(as.double(temp_t_result$statistic))
  }

  ########################################### 
  #  Calculate differential correlation
  ########################################### 
  print("Calculating pairwise correlation for normal states...")
  # get correlation for normal group
  r_normal =cor(t(expressionMatrix_normal),use ="pairwise.complete.obs")
  # increase memory size
  # memory.limit(10000)
  r_normal[which(r_normal ==1)]=0.9999999  
  # get differential co-expression, z, for normal group
  gc()  
  z_1=(1+r_normal)
  gc()  
  z_1=z_1/(1-r_normal)
  rm(r_normal)
  gc()  
  z_1=abs(z_1)
  gc() 
  z_1=0.5*log(z_1)
  gc()  

  print("Calculating pairwise correlation for disease states...")
  # get correlation for disease group
  r_disease =cor(t(expressionMatrix_disease),use ="pairwise.complete.obs")
  r_disease[which(r_disease ==1)]=0.9999999
  # get differential co-expression, z, for disease group
  gc()  
  z_2=(1+r_disease)
  gc()  
  z_2=z_2/(1-r_disease)
  rm(r_disease)
  gc()  
  z_2=abs(z_2)
  gc() 
  z_2=0.5*log(z_2)
  gc()  

  print("Calculating differential co-expression measures ...")
  n1=ncol(expressionMatrix_normal)
  n2=ncol(expressionMatrix_disease) 
  d_r =abs(z_1-z_2)/sqrt(1/(n1-3)+1/(n2-3)) 

  # save memory
  rm(z_1)
  rm(z_2)
  gc() 
  
  ########################################### 
  # get pathway data
  ########################################### 
  minSupport =3
  print("Reading functional gene set data")
  pathway = getPathway(geneSetInputFile,geneName,minSupport)  

  ########################################### 
  # get optimal threshold by maximizing chi-square value
  ########################################### 
  # import file for calculating FDR
  # import file for calculating adjusted p-values
  print("Identifying optimal thresholds for genes")
  optimalCutOff= getDE_DC_OptimalThreshold(t_result, MaxGene, d_r,minSupport)
  #save.image("image_liver")

  ########################################### 
  # identify best associated functional gene sets for each partitions
  ########################################### 
  outputFileName = "temp"
 
  bestAssoGeneSet_HC_HD=data.frame(index= double(),name=character(), obsA= double(), expA= double(), pValue= double(), matchedGene=list(), stringsAsFactors = FALSE,LD= double(), HD= double(), LC= double(), HC= double(), LC_LD= double(), LC_HD= double(), HC_LD= double(), HC_HD= double())

  # import files to identify best associated functional gene sets for each partitions
  print("Identifying best associated functional gene set for each gene...")
  # for every gene i
  for (i in 1:MaxGene) {
    if (i %% 1==0) {
      print(sprintf("Gene id: %d",i))
    }
    # identify genes in each of the partitions
    one_8partitions = getPartitionIndex(i,t_result,optimalCutOff,d_r)
    
  # identify the best associated gene sets for HDC-HDE partitions
    tempBestAssoGeneSet = getBestAssociatedGeneSet(pathway, one_8partitions, one_8partitions$highR_highDE[[1]],MaxGene,minSupport)
    bestAssoGeneSet_HC_HD =rbind(bestAssoGeneSet_HC_HD,tempBestAssoGeneSet)

  }  

###########################################   
# perform p-value adjustments (FDR) for each p-value columns
# input: best associated gene set with unadjusted p-values
# output: best associated gene set with adjusted p-values
########################################### 
getFDRFor9Columns = function(tempAssoGeneSet) {
  tempAssoGeneSet[,"pValue"] = getFDR(tempAssoGeneSet[,"pValue"])
  tempAssoGeneSet[,"FI_LD"] = -1*log(getFDR(tempAssoGeneSet[,"FI_LD"]),2)
  tempAssoGeneSet[,"FI_HD"] = -1*log(getFDR(tempAssoGeneSet[,"FI_HD"]),2)
  tempAssoGeneSet[,"FI_LC"] = -1*log(getFDR(tempAssoGeneSet[,"FI_LC"]),2)
  tempAssoGeneSet[,"FI_HC"] = -1*log(getFDR(tempAssoGeneSet[,"FI_HC"]),2)
  tempAssoGeneSet[,"FI_LC_LD"] = -1*log(getFDR(tempAssoGeneSet[,"FI_LC_LD"]),2)
  tempAssoGeneSet[,"FI_LC_HD"] = -1*log(getFDR(tempAssoGeneSet[,"FI_LC_HD"]),2)
  tempAssoGeneSet[,"FI_HC_LD"] = -1*log(getFDR(tempAssoGeneSet[,"FI_HC_LD"]),2)
  tempAssoGeneSet[,"FI_HC_HD"] = -1*log(getFDR(tempAssoGeneSet[,"FI_HC_HD"]),2)
  return (tempAssoGeneSet)
}


  # adjust p-values in multiple tests
  bestAssoGeneSet_HC_HD = getFDRFor9Columns(bestAssoGeneSet_HC_HD)

########################################### 
# function to print one row of the results for the gene i's partitions 
# input: (1) significant level for best associated gene set
#        (2) best associated gene sets for the partitions of gene i
# output: result to text file
########################################### 
printOnePartitionResult = function(significanceLevel,tempAssoGeneSet) {
  tempPValue= tempAssoGeneSet[1,"pValue"]
  
  # if gene set is significant 
  if (tempPValue<significanceLevel) {
    tempGeneSetIndex = tempAssoGeneSet[1,"index"]


    cat(sprintf("%s\t",pathway[tempGeneSetIndex,"name"]))
    cat(sprintf("%s\t",pathway[tempGeneSetIndex,"categoryInfo"]))
    cat(sprintf("%d\t",tempAssoGeneSet[1,"obsA"]))
    cat(sprintf("%.7f\t",tempAssoGeneSet[1,"expA"]))
    cat(sprintf("%e\t",tempPValue))
    tempMatchedGenesIndex = tempAssoGeneSet$matchedGene[[1]]
    tempMatchedGenesIndexCount =length(tempMatchedGenesIndex)
    for (j in 1:tempMatchedGenesIndexCount) {
      cat(sprintf("%s ",geneName[tempMatchedGenesIndex[j]]))
    }
    cat(sprintf("\t"))
    cat(sprintf("%.7f\t",tempAssoGeneSet[1,"FI_LD"]))
    cat(sprintf("%.7f\t",tempAssoGeneSet[1,"FI_HD"]))
    cat(sprintf("%.7f\t",tempAssoGeneSet[1,"FI_LC"]))
    cat(sprintf("%.7f\t",tempAssoGeneSet[1,"FI_HC"]))
    cat(sprintf("%.7f\t",tempAssoGeneSet[1,"FI_LC_LD"]))
    cat(sprintf("%.7f\t",tempAssoGeneSet[1,"FI_LC_HD"]))
    cat(sprintf("%.7f\t",tempAssoGeneSet[1,"FI_HC_LD"]))
    cat(sprintf("%.7f\t",tempAssoGeneSet[1,"FI_HC_HD"]))
  } else {
  # if gene set is not significant 
    cat("-\t")
    cat("-\t")
    cat("0\t")
    cat("0\t")
    cat("1\t")
    cat("-\t")
    cat("0\t")
    cat("0\t")
    cat("0\t")
    cat("0\t")
    cat("0\t")
    cat("0\t")
    cat("0\t")
    cat("0\t")
  }
}
  ########################################### 
  # open file and write results
  ########################################### 
  openFileToWrite(outputFileName)
  print("Processing raw results...")
  
  sink(outputFileName,append=TRUE)
  for (i in 1:MaxGene) {
    cat(sprintf("%d\t",i))
    cat(sprintf("%s\t",geneName[i]))
    cat(sprintf("%.7f\t",t_result[i,"tScore"]))
    cat(sprintf("%.7f\t",t_result[i,"absTScore"]))
    cat(sprintf("%e\t",t_result[i,"pValue"]))
    cat(sprintf("%.7f\t",foldChange[i]))

    cat(sprintf("%.7f\t",optimalCutOff[i,"tCutOff"]))
    cat(sprintf("%.7f\t",optimalCutOff[i,"rCutOff"]))
    cat(sprintf("%d\t",optimalCutOff[i,"obsA"]))
    cat(sprintf("%d\t",optimalCutOff[i,"obsB"]))
    cat(sprintf("%d\t",optimalCutOff[i,"obsC"]))
    cat(sprintf("%d\t",optimalCutOff[i,"obsD"]))
    cat(sprintf("%.7f\t",optimalCutOff[i,"expA"]))
    cat(sprintf("%.7f\t",optimalCutOff[i,"chi2Value"]))
    cat(sprintf("%e\t",optimalCutOff[i,"pValue"]))
    
    printOnePartitionResult(significanceLevel,bestAssoGeneSet_HC_HD[i,])
    cat("\n")  
  }
  sink()  

  sumResult_MinGain()
}
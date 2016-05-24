########################################################################## 
#' Summarize the functional gene set results into text file
#' @export
########################################################################## 
# DECODE (c) Copyright 2014 by The Hong Kong Polytechnic University, Department of Health Technology and Informatics 
# Written by Thomas Lui
# Permission is granted to copy and use this program provided no fee is
# charged for it and provided that this copyright notice is not removed.
##########################################################################
sumResult_MinGain = function() {
  rm(list = ls())
  numberOfTopGeneSelected=5

  ###########################
  # read input data 
  ###########################
  print("Summarizing functional gene set results...")
  rawData = read.table("temp", sep="\t",row.names = NULL, quote = "", header=TRUE,fill =FALSE)  
	
  outputFileName = "out_summary.txt"
  threshold_optimalChi=0.05
  minSupport=1
  

  dim_rawData=dim(rawData)
  MaxRow =nrow(rawData)
  MaxCol =ncol(rawData)
  optimalChi =data.frame(geneName=character(), obsA=double(), expA=double(), chi2Value=double(), chiPValue=double(), stringsAsFactors = FALSE)
  optimalChi [1:MaxRow,"geneName"]=as.matrix(rawData[1:MaxRow,"Gene"])
  ##############
  # normal
  ############## 
  tempCol=9
  optimalChi [1:MaxRow,"obsA"]=as.matrix(rawData[1:MaxRow,tempCol])   
  optimalChi [1:MaxRow,"expA"]=as.matrix(rawData[1:MaxRow,tempCol+4]) 
  optimalChi [1:MaxRow,"chi2Value"]=as.matrix(rawData[1:MaxRow,tempCol+5])
  optimalChi [1:MaxRow,"chiPValue"]=as.matrix(rawData[1:MaxRow,tempCol+6])
  
  highR_highDE_partition = data.frame(geneSet=character(), geneSetCategory=character(), pValue=double(), FI_LD=double(), FI_HD=double(), FI_LR=double(),FI_HR=double(),FI_LR_LD=double(),FI_LR_HD=double(),FI_HR_LD=double(),FI_HR_HD=double(), stringsAsFactors = FALSE)
  tempCol=tempCol+7
  highR_highDE_partition[1:MaxRow,"geneSet"]=as.matrix(rawData[1:MaxRow,tempCol])
  highR_highDE_partition[1:MaxRow,"geneSetCategory"]=as.matrix(rawData[1:MaxRow,tempCol+1])
  highR_highDE_partition[1:MaxRow,"pValue"]=as.matrix(rawData[1:MaxRow,tempCol+4])
  highR_highDE_partition[1:MaxRow,"FI_LD"]=as.matrix(rawData[1:MaxRow,tempCol+6])
  highR_highDE_partition[1:MaxRow,"FI_HD"]=as.matrix(rawData[1:MaxRow,tempCol+7])
  highR_highDE_partition[1:MaxRow,"FI_LR"]=as.matrix(rawData[1:MaxRow,tempCol+8])
  highR_highDE_partition[1:MaxRow,"FI_HR"]=as.matrix(rawData[1:MaxRow,tempCol+9])
  highR_highDE_partition[1:MaxRow,"FI_LR_LD"]=as.matrix(rawData[1:MaxRow,tempCol+10])
  highR_highDE_partition[1:MaxRow,"FI_LR_HD"]=as.matrix(rawData[1:MaxRow,tempCol+11])
  highR_highDE_partition[1:MaxRow,"FI_HR_LD"]=as.matrix(rawData[1:MaxRow,tempCol+12])
  highR_highDE_partition[1:MaxRow,"FI_HR_HD"]=as.matrix(rawData[1:MaxRow,tempCol+13])
  highR_highDE_partition[1:MaxRow,"Gain_over_DE"]=(as.matrix(rawData[1:MaxRow,tempCol+13]) - as.matrix(rawData[1:MaxRow,tempCol+7])) 
  highR_highDE_partition[1:MaxRow,"Gain_over_DC"]=(as.matrix(rawData[1:MaxRow,tempCol+13])- as.matrix(rawData[1:MaxRow,tempCol+9]))
  highR_highDE_partition[1:MaxRow,"Min_Gain"]=apply(cbind(highR_highDE_partition[1:MaxRow,"Gain_over_DE"], highR_highDE_partition[1:MaxRow,"Gain_over_DC"]),1,min)
  #unlink("temp")

###########
# summary function for each partition
############  
getSummary = function(optimalChi, partition, signOfAsso,outputFileName,tableDescription,rankedGeneSet) {
  write(tableDescription,outputFileName,append=TRUE)
  partition= cbind(optimalChi,partition)
  if (signOfAsso==0){
    tempIndex = partition$chiPValue<threshold_optimalChi & partition$pValue<threshold_optimalChi 
  }
  if (signOfAsso==1) {
    tempIndex = partition$obsA > partition$expA &  partition$chiPValue<threshold_optimalChi & partition$pValue<threshold_optimalChi 
  } 
  if (signOfAsso==-1) {
    tempIndex = partition$obsA < partition$expA &  partition$chiPValue<threshold_optimalChi  & partition$pValue<threshold_optimalChi
  }  
  colCount =ncol(partition)
  selectedPartition = partition[tempIndex,1:colCount]

  uniqueGeneSet = unique(selectedPartition[,"geneSet"]) # get the gene set name
  uniqueGeneSet = uniqueGeneSet[which(!is.na(uniqueGeneSet))]  # remove NULL
  uniqueGeneSetCount= length(uniqueGeneSet)
  
  summary_selectedPartition = data.frame(geneSet=character(),geneSetCategory=character(),geneSetCount=double(),meanFI_for_HDE=double(), meanFI_for_HDC=double(),meanFI_for_HDC_HDE=double(), FIGain_over_DE =double(), FIGain_over_DC =double(), MinFI_Gain=double(), stringsAsFactors = FALSE)
  if (uniqueGeneSetCount>0) {
    for (i in 1:uniqueGeneSetCount) {

      summary_selectedPartition[i,"geneSet"]=uniqueGeneSet[i]
	  tempIndex =which(selectedPartition$geneSet == uniqueGeneSet[i])
      summary_selectedPartition[i,"geneSetCategory"]=selectedPartition$geneSetCategory[tempIndex[1]]
	  summary_selectedPartition[i,"geneSetCount"] =length(tempIndex)
      summary_selectedPartition[i,"meanFI_for_HDE"]= mean(selectedPartition$FI_HD[tempIndex])
      summary_selectedPartition[i,"meanFI_for_HDC"]= mean(selectedPartition$FI_HR[tempIndex])
      summary_selectedPartition[i,"meanFI_for_HDC_HDE"]= mean(selectedPartition$FI_HR_HD[tempIndex])
      summary_selectedPartition[i,"FIGain_over_DE"]= mean(selectedPartition$Gain_over_DE[tempIndex])
      summary_selectedPartition[i,"FIGain_over_DC"]= mean(selectedPartition$Gain_over_DC[tempIndex])
      summary_selectedPartition[i,"MinFI_Gain"]= mean(selectedPartition$Min_Gain[tempIndex])
	 
	  
      #rank by Min_Gain
      tempRankedGeneIndex =rank(selectedPartition$Min_Gain[tempIndex], ties.method = "first")
      tempRankedGene=0
      tempRankedGene[tempRankedGeneIndex] =selectedPartition$geneName[tempIndex]
      tempRankedGene=rev(tempRankedGene)
      tempRankedGeneCount=length(tempRankedGene)  

	}

    #rank by min FI gain
    summary_selectedPartition=summary_selectedPartition[order(summary_selectedPartition$MinFI_Gain, decreasing=TRUE),]
	if (length(rankedGeneSet)==0) {
	} else {
	  # for unmatched gene sets sort according to number of gene set count
	  uniqueGeneSet= summary_selectedPartition[,"geneSet"]	  
	  newOrderIndex = match(rankedGeneSet,uniqueGeneSet, nomatch=-1)
      newOrderIndex=newOrderIndex[newOrderIndex!=-1]
      uniqueGeneSetNotInRankedGeneSet = match(uniqueGeneSet,rankedGeneSet, nomatch=-1)
	  uniqueGeneSetNotInRankedGeneSet =uniqueGeneSet[uniqueGeneSetNotInRankedGeneSet ==-1]
	  uniqueGeneSetNotInRankedGeneSetIndex = match(uniqueGeneSetNotInRankedGeneSet,uniqueGeneSet, nomatch=-1)
	  newOrderIndex =c(newOrderIndex,uniqueGeneSetNotInRankedGeneSetIndex)
      summary_selectedPartition=summary_selectedPartition[newOrderIndex,]
	}


    summary_selectedPartition =summary_selectedPartition[summary_selectedPartition$geneSetCount>=minSupport,]
    suppressWarnings(write.table(summary_selectedPartition, outputFileName, sep='\t', quote=FALSE,row.names = FALSE, append=TRUE))
  }
  if (length(rankedGeneSet)==0) {
    return(summary_selectedPartition[,"geneSet"])
  }
}

  write("",outputFileName,append=FALSE)
  
  rankedGeneSet=vector()
  rankedGeneSet = getSummary(optimalChi,highR_highDE_partition,1,outputFileName, "\nBest associated gene sets with highest mean minimum functional information (FI) gain for HDC_HDE partitions\n",rankedGeneSet)
  print("Done. Result is saved in out_summary.txt")
}
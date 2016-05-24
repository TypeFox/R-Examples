#####################################################
#' Get best associated functional gene sets for partitions of gene i
#' @param pathway All functional gene sets
#' @param all8Partitions All eight possible partitions for gene i
#' @param onePartition The partition to be associated with the functional gene set
#' @param MaxGene Number of genes in expression data
#' @param minSupport Minimum support for functional gene set
#' @return The adjusted p-values for the best associated gene set of the input partition
#' @export
#####################################################
# DECODE (c) Copyright 2014 by The Hong Kong Polytechnic University, Department of Health Technology and Informatics 
# Written by Thomas Lui
# Permission is granted to copy and use this program provided no fee is
# charged for it and provided that this copyright notice is not removed.
##########################################################################
getBestAssociatedGeneSet = function(pathway, all8Partitions, onePartition ,MaxGene,minSupport) {
  optimal_tempAGenesIndex= vector()
  onePartitionCount = length(onePartition)
  pathwayCount = nrow(pathway)
  #########################################
  #  Find max associated gene set
  #########################################
  ########################################################
  # association with pathway sets
  #                    in gene set ,  not in gene set
  # gene list             A               B
  # not in gene list      C               D
  #######################################################
  optimal_pathwayIndex=0
  optimal_tempA=0
  optimal_tempB=0
  optimal_tempC=0
  optimal_tempD=0
  optimal_expectedA=0
  optimal_pValue=1

  tempPathIndexSetSize =pathway[,"size_indexSet"]
  multipleTestCount=0
  for (j in 1:pathwayCount) {
    ##############################################
    # get overlap of 2 lists
    ##############################################
	tempAGenesIndex =intersect(pathway$indexSet[[j]],onePartition) 
    tempA =length(tempAGenesIndex)
	if (tempA < minSupport) {
	  next
	}
    tempB = onePartitionCount - tempA
    tempC = tempPathIndexSetSize[j] - tempA
    expectedA = (tempA+tempB)*(tempA+tempC)/ MaxGene
	if (tempA<expectedA) {
	  next
	}
    tempD = MaxGene-onePartitionCount - tempC
    tempConTable = matrix(c(tempA, tempC , tempB, tempD),2,2)
    ##############################################
    #get p-values
    ##############################################
    fisher_result = fisher.test(tempConTable)
    pValue = fisher_result$p.value
	multipleTestCount= multipleTestCount+1
    if (pValue<optimal_pValue) {
      optimal_pathwayIndex=j
	  optimal_pValue=pValue
      optimal_tempA=tempA
      optimal_tempB=tempB
      optimal_tempC=tempC
      optimal_tempD=tempD
      optimal_expectedA=expectedA
	  optimal_tempAGenesIndex =tempAGenesIndex
 	}
  }
  optimal_pValue = optimal_pValue * multipleTestCount # Bonferron correction
  if (optimal_pValue>1) {
    optimal_pValue=1
	optimal_pathwayIndex=0
  }
  if (optimal_pValue<10^-309) {   # set it to zero as excel (64-bits) can't handle it anyway
    optimal_pValue=0
  }
  bestAssoGeneSet=data.frame(index= double(),name=character(), obsA= double(), expA= double(), pValue= double(), matchedGene=list(), stringsAsFactors = FALSE,
                             FI_LD= double(), FI_HD= double(), FI_LC= double(), FI_HC= double(), FI_LC_LD= double(), FI_LC_HD= double(), FI_HC_LD= double(), FI_HC_HD= double())
  if (optimal_pathwayIndex !=0) {
	bestAssoGeneSet[1,"index"] =optimal_pathwayIndex
	bestAssoGeneSet[1,"name"] =pathway[optimal_pathwayIndex,"name"]
	bestAssoGeneSet[1,"obsA"] =optimal_tempA
	bestAssoGeneSet[1,"expA"] =optimal_expectedA
	bestAssoGeneSet[1,"pValue"] =optimal_pValue
	bestAssoGeneSet$matchedGene[[1]] =optimal_tempAGenesIndex
    
    # find the p-value for this gene set in other partitions
    bestAssoGeneSet[1,"FI_LD"] = getAssoGeneSetPValue(all8Partitions$lowDE[[1]], pathway$indexSet[[optimal_pathwayIndex]],multipleTestCount,MaxGene)
    bestAssoGeneSet[1,"FI_HD"] = getAssoGeneSetPValue(all8Partitions$highDE[[1]], pathway$indexSet[[optimal_pathwayIndex]],multipleTestCount,MaxGene)
    bestAssoGeneSet[1,"FI_LC"] = getAssoGeneSetPValue(all8Partitions$lowR[[1]], pathway$indexSet[[optimal_pathwayIndex]],multipleTestCount,MaxGene)
    bestAssoGeneSet[1,"FI_HC"] = getAssoGeneSetPValue(all8Partitions$highR[[1]], pathway$indexSet[[optimal_pathwayIndex]],multipleTestCount,MaxGene)
    bestAssoGeneSet[1,"FI_LC_LD"] = getAssoGeneSetPValue(all8Partitions$lowR_lowDE[[1]], pathway$indexSet[[optimal_pathwayIndex]],multipleTestCount,MaxGene)
    bestAssoGeneSet[1,"FI_LC_HD"] = getAssoGeneSetPValue(all8Partitions$lowR_highDE[[1]], pathway$indexSet[[optimal_pathwayIndex]],multipleTestCount,MaxGene)
    bestAssoGeneSet[1,"FI_HC_LD"] = getAssoGeneSetPValue(all8Partitions$highR_lowDE[[1]], pathway$indexSet[[optimal_pathwayIndex]],multipleTestCount,MaxGene)
    bestAssoGeneSet[1,"FI_HC_HD"] = getAssoGeneSetPValue(all8Partitions$highR_highDE[[1]], pathway$indexSet[[optimal_pathwayIndex]],multipleTestCount,MaxGene)

  } else {
	bestAssoGeneSet[1,"index"] =-1
	bestAssoGeneSet[1,"name"] ="NA"
	bestAssoGeneSet[1,"obsA"] =0
	bestAssoGeneSet[1,"expA"] =0
	bestAssoGeneSet[1,"pValue"] =1
	bestAssoGeneSet$matchedGene[[1]] =-1
    
    # find the p-value for this gene set in other partitions	
    bestAssoGeneSet[1,"FI_LD"] = 1
    bestAssoGeneSet[1,"FI_HD"] = 1
    bestAssoGeneSet[1,"FI_LC"] = 1
    bestAssoGeneSet[1,"FI_HC"] = 1
    bestAssoGeneSet[1,"FI_LC_LD"] = 1
    bestAssoGeneSet[1,"FI_LC_HD"] = 1
    bestAssoGeneSet[1,"FI_HC_LD"] = 1
    bestAssoGeneSet[1,"FI_HC_HD"] = 1
  }
  return(bestAssoGeneSet)

}
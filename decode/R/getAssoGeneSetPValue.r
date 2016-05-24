########################################################################## 
#' Calculate the p-value between selected genes and functional gene set
#' @param geneList Selected genes 
#' @param geneSet Functional gene set
#' @param multipleTestCount Number of multiple testing
#' @param MaxGene Number of genes in expression data
#' @return The adjusted p-value for the associated gene set
#' @export
########################################################################### 
# DECODE  (c) Copyright 2014 by The Hong Kong Polytechnic University, Department of Health Technology and Informatics 
# Written by Thomas Lui
# Permission is granted to copy and use this program provided no fee is
# charged for it and provided that this copyright notice is not removed.
##########################################################################
getAssoGeneSetPValue= function(geneList ,geneSet, multipleTestCount,MaxGene) {
  geneListCount = length(geneList)
  geneSetCount= length(geneSet)
  ########################################################
  # association with pathway sets
  #                    in gene set ,  not in gene set
  # gene list             A               B
  # not in gene list      C               D
  #######################################################

  # get overlap of 2 lists
  tempAGenesIndex =intersect(geneList,geneSet)
  tempA =length(tempAGenesIndex)
  tempB = geneListCount - tempA
  tempC = geneSetCount - tempA
  expectedA = (tempA+tempB)*(tempA+tempC)/ MaxGene
  if (tempA<expectedA) {  # consider only over representation
    return (1)
  }
  tempD = MaxGene-geneListCount - tempC
  tempConTable = matrix(c(tempA, tempC , tempB, tempD),2,2)

  ##############################################
  #get p-values
  ##############################################
  fisher_result = fisher.test(tempConTable)
  pValue = fisher_result$p.value
  pValue = pValue* multipleTestCount  # Bonferron correction
  if (pValue>1) {
    pValue =1
  }
  if (pValue<10^-309) {   # set it to zero as excel (64-bits) can't handle it anyway
    pValue =0
  }
  
  return(pValue) 
}
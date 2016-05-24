############################################################
#' read functional gene sets
#' @param inputFile Input file name
#' @param geneName Gene name lists
#' @param minSupport Minimum support
#' @return Functional gene set
#' @export
############################################################
# DECODE  (c) Copyright 2014 by The Hong Kong Polytechnic University, Department of Health Technology and Informatics 
# Written by Thomas Lui
# Permission is granted to copy and use this program provided no fee is
# charged for it and provided that this copyright notice is not removed.
##########################################################################
getPathway = function(inputFile,geneName,minSupport) {
  inFile  = file(inputFile, open = "r")
  pathwayName="";
  pathwayCategory="";
  pathwayIndexSet=list();
  pathwayCount=0;
  size_pathwayIndexSet=0;

  pathway=data.frame(name=character(), categoryInfo=character(),indexSet=list(), size_indexSet=double(),stringsAsFactors = FALSE)
  
  
  while (TRUE) {
    # read a line
    oneLine = readLines(inFile, n = 1, warn = FALSE) 
	if (length(oneLine)== 0) {
	  break
	}
	if (pathwayCount %% 100 ==0) {
      # print(pathwayCount)
    }
    oneLineArray = unlist(strsplit(oneLine, "\t"))

    onePathwayName =  oneLineArray[1]
    onePathwayCategory=  oneLineArray[2]
	size_oneLineArray = length(oneLineArray)
	
	tempPathwayIndexSetSize=0
    onePathwayGeneSet=oneLineArray[3:size_oneLineArray]

	onePathwayGeneSetIndex = match(onePathwayGeneSet,geneName, nomatch=-1)
	onePathwayGeneSetIndex =onePathwayGeneSetIndex[onePathwayGeneSetIndex != -1]
	size_onePathwayGeneSetIndex = length(onePathwayGeneSetIndex)
	#at least there are some matches (>=minSupport) of this pathway set to the genes of the gene expression data
	if (size_onePathwayGeneSetIndex >=minSupport) {
      pathwayCount=pathwayCount+1
      pathway[pathwayCount,"name"]= onePathwayName
      pathway[pathwayCount,"categoryInfo"]=onePathwayCategory
	  pathway$indexSet[[pathwayCount]] =onePathwayGeneSetIndex  # assign list of numbers
	  pathway[pathwayCount,"size_indexSet"]=size_onePathwayGeneSetIndex
	}
  } 
  close(inFile)

  return(pathway)
}
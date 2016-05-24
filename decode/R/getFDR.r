#####################################################
#' Adjust p-value by Benjamini and Hochberg method
#' @param pValues Unadjusted p-values
#' @return Adjusted p-values
#' @export
#####################################################
# DECODE (c) Copyright 2014 by The Hong Kong Polytechnic University, Department of Health Technology and Informatics 
# Written by Thomas Lui
# Permission is granted to copy and use this program provided no fee is
# charged for it and provided that this copyright notice is not removed.
##########################################################################
getFDR = function(pValues) {
  tempRank = rank(pValues, ties.method= "min")
  numOfPvalue = length(pValues)
  adjustedPValues =pValues*numOfPvalue/tempRank
  adjustedPValues[which(adjustedPValues>1)]=1  # equal 1 if adjusted value greater than 1
  return(adjustedPValues)
}
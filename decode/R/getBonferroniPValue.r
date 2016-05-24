#####################################################
#' Adjust p-value by Bonferroni correction
#' @param pValues Unadjusted p-values
#' @return Adjusted p-values
#' @export
#####################################################
# DECODE (c) Copyright 2014 by The Hong Kong Polytechnic University, Department of Health Technology and Informatics 
# Written by Thomas Lui
# Permission is granted to copy and use this program provided no fee is
# charged for it and provided that this copyright notice is not removed.
##########################################################################
getBonferroniPValue = function(pValues) {
  numOfPvalue = length(pValues)
  adjustedPValues =pValues*numOfPvalue
  adjustedPValues[which(adjustedPValues>1)]=1  # equal 1 if adjusted value greater than 1
  return(adjustedPValues)
}
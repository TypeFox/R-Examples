############################################
#' Get gene index of 8 partitions for gene i
#' @param gene_i Gene i index
#' @param t_result t-statistics
#' @param optimalCutOff Optimal thresholds
#' @param abs_r Matrix consisting of absolute values of all differential co-expression measures
#' @return The selected genes for each partition in index
#' @export
################################################
# DECODE v1.0  (c) Copyright 2014 by The Hong Kong Polytechnic University, Department of Health Technology and Informatics 
# Written by Thomas Lui
# Permission is granted to copy and use this program provided no fee is
# charged for it and provided that this copyright notice is not removed.
##########################################################################
getPartitionIndex = function(gene_i,t_result, optimalCutOff, abs_r) {
	
  # select genes greater or equal to the optimal t
  clusterGeneIndex_highDE = t_result$absTScore>=optimalCutOff[gene_i,"tCutOff"]
  # select genes lower the optimal t
  clusterGeneIndex_lowDE = !clusterGeneIndex_highDE	  
  # select genes greater or equal to the optimal dc
  clusterGeneIndex_highR = abs_r[gene_i,]>=optimalCutOff[gene_i,"rCutOff"]
  # select genes lower the optimal dc
  clusterGeneIndex_lowR = !clusterGeneIndex_highR
  
  clusterGeneIndex_lowR_lowDE = clusterGeneIndex_lowR & clusterGeneIndex_lowDE
  clusterGeneIndex_lowR_highDE = clusterGeneIndex_lowR & clusterGeneIndex_highDE
  clusterGeneIndex_highR_lowDE =clusterGeneIndex_highR & clusterGeneIndex_lowDE
  clusterGeneIndex_highR_highDE = clusterGeneIndex_highR & clusterGeneIndex_highDE
  
  partition=data.frame(id= double(),lowDE=list(), highDE=list(), lowR=list(), highR=list(), lowR_lowDE=list(), lowR_highDE=list(), highR_lowDE=list(), highR_highDE=list(), stringsAsFactors = FALSE)
  partition[1,"id"]=1
  partition$lowDE[[1]] = which(clusterGeneIndex_lowDE)
  partition$highDE[[1]] = which(clusterGeneIndex_highDE)
  partition$lowR[[1]] = which(clusterGeneIndex_lowR)
  partition$highR[[1]] = which(clusterGeneIndex_highR)
  partition$lowR_lowDE[[1]] = which(clusterGeneIndex_lowR_lowDE)
  partition$lowR_highDE[[1]] = which(clusterGeneIndex_lowR_highDE)
  partition$highR_lowDE[[1]] = which(clusterGeneIndex_highR_lowDE)
  partition$highR_highDE[[1]] = which(clusterGeneIndex_highR_highDE)
	  
  return(partition)
}
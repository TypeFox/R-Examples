#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


shockSelect <- function(expdata){
  
  if(is.matrix(expdata) == FALSE & is.data.frame(expdata) == FALSE) 
    stop(paste(sQuote("expdata"), "must be a matrix"))
  
  ## function shockSelect
  ## input expdata (matrix): data 
  ## output SHDJlabels (vector): vector of partition labels based on the slope heuristic dimension jump
  ## output SHDJlabels (vector): vector of partition labels based on the slope heuristic robust regression
  ## output capusheOuput (capushe objet): output of the $\kappa$ coefficient calibration function capushe 
  
  ## require
  ## library(mvtnorm)
  ## library(capushe)
  ## library(igraph)
  
  
  ## dimension of the dataset expdata
  n <- dim(expdata)[1]
  p <- dim(expdata)[2]
  
  
  ## threshold of absS matrix
  myPartition <- thresholdAbsSPath(expdata)$partitionList
  maxCC <- unlist(lapply(myPartition, function(x) max(table(x))))
  resLogLike <- lapply(myPartition, function(x) computeLoglikeFromPartition(x,expdata))
  myLog <- unlist(lapply(resLogLike, function(x) x$loglike))
  myDf <- unlist(lapply(resLogLike, function(x) x$df))
  
  
  ## SLOPE HEURISTIQUE SELECTION
  ## first, we remove points with infinite log-likelihood
  myDf.clean <- myDf[myLog!="-Inf"]
  myLog.clean <- myLog[myLog!="-Inf"]
  myPartition.clean <- myPartition[myLog!="-Inf"]
  
  ## we remove points with 0 degree of freedom 
  mat <- cbind(as.vector(myDf.clean[myDf.clean>0]), as.vector(myDf.clean[myDf.clean>0]), as.vector(myDf.clean[myDf.clean>0]), -as.vector(myLog.clean[myDf.clean>0]))
  ResCapushe <- capushe(mat, n, psi.rlm="lm", pct=0.10)
  ## we collect the labels for each selected partition
  SHDJlabels <- myPartition.clean[[min(which(myDf.clean==as.numeric(ResCapushe@Djump@model)))]]
  SHRRlabels<- myPartition.clean[[min(which(myDf.clean==as.numeric(ResCapushe@DDSE@model)))]]
  
  return(list(SHDJlabels = SHDJlabels, SHRRlabels = SHRRlabels, capusheOutput = ResCapushe ))
}




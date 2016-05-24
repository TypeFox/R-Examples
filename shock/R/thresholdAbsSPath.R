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




thresholdAbsSPath <- function(expdata){
  
  ## function that take data
  ## and provide the clustering output
  ## based on the threshold sample covariance matrix
  ## input expdata (matrix): data to compute the sample covariance matrix S
  ## output partitionList (list): list of resulting partitions from thresolding S
  ## output lambdaPath (list): list of threshold parameters
  
  ## require
  ## library(igraph)
  
  
  if(is.matrix(expdata) == FALSE & is.data.frame(expdata) == FALSE) 
    stop(paste(sQuote("expdata"), "must be a matrix"))
   
  labelsPath <- list()
  Sabs <- abs(cor(expdata))
  valThres <- Sabs[upper.tri(Sabs)]
  orderValue <- valThres[order(valThres)]
  
  for( lam in 1:length(orderValue)){
    lambdaR <- orderValue[lam]
    E <- Sabs
    E[Sabs>lambdaR] <- 1
    E[Sabs<lambdaR] <- 0
    E[Sabs==lambdaR] <- 0
    goutput <- graph.adjacency(E,mode="undirected",weighted=NULL)
    labelsPath[[lam]] <- clusters(goutput)$membership    
  }
  
  return(list(partitionList=unique(labelsPath), lambdaPath=unique(orderValue)))
}



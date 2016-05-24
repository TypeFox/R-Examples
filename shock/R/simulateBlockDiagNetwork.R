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




simulateBlockDiagNetwork <- function (p, labels)
{
  
  if(is.list(labels) == FALSE & is.vector(labels) == FALSE & is.factor(labels) == FALSE) 
    stop(paste(sQuote("labels"), "must be a list"))
  
  
  ## function to simulate Graph
  ## modif from GGMselect simulateNetwork
  ## require GGMselect for simulateGraph
  
  ## input p (integer): number of variables input
  ## input labels (factor vector): labels index for each p variables
  ## input eta (real) 
  ## input extraeta (real)
  ## eta level of connectedness inside modules from GGMselect:::simulateGraph
  ## extraeta level of connectedness outside modules from GGMselect:::simulateGraph
  ## to simulate the network, simulateGraph use split the p variables into 3 sets
  ## put edges between subsets with probability extraeta
  ## put edges within subsets with probability eta
  ## eta and extraeta are for each block
  
  
  
  ## adjacency matrix
  A <- matrix(0, p, p)
  ## covariance matrix
  C <- matrix(0,p,p)
  ## partial correlation matrix
  PCor <- matrix(0,p,p)
  ## number of labels
  Q <- length(unique(labels))
  
  
  ## sample edges
  for (q in 1:Q) {
    ps <- length(labels[labels == q])
    
    graphSmall <- simulateGraph(ps,eta=1,extraeta=1)
    ## code to test if block are really connected component
    ## g <- graph.adjacency(graphSmall$G,mode="undirected",weighted=NULL)
    ## if(labels(g)$no !=1) warning("Connected component split")
    
    A[labels == q, labels == q] <- graphSmall$G
    C[labels == q, labels == q] <- graphSmall$C
    PCor[labels == q, labels == q] <- graphSmall$PCor
    
  }
   
  nodes <- as.character(paste("g", 1:p, sep = ""))
  dimnames(A) <- list(nodes, nodes)
  dimnames(C) <- list(nodes, nodes)
  dimnames(PCor) <- list(nodes, nodes)
  return(structure(list(A = A, C = C, PCor = PCor, labels = labels), class = "block.diag.network"))
}



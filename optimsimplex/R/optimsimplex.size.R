# Copyright (C) 2008-2009 - INRIA - Michael Baudin
# Copyright (C) 2009-2010 - DIGITEO - Michael Baudin
# Copyright (C) 2010-2014 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#
# This source code is a R port of the optimsimplex component
# originally written by Michael Baudin for Scilab.

optimsimplex.size <- function(this=NULL,method=NULL){

  n <- this$n
  nv <- this$nbve
  if (is.null(method)) method <- 'sigmaplus'

  if (!any(method==c('Nash','sigmaplus','sigmaminus')))
    stop(sprintf('optimsimplex.size: Unknown simplex size method %s', method))

  if (method=='Nash'){
    v1 <- matrix(rep(this$x[1,],nv-1),nrow=nv-1,byrow=TRUE)
    edges <- this$x[2:nv,,drop=FALSE] - v1
    abedges <- abs(edges)
    n1 <- apply(abedges,1,sum)
    ssize <- sum (n1)
  }
  if (method=='diameter'){
    ssize <- 0.0
    for (i in 1:nv){
      vi <- matrix(rep(this$x[i,],nv),nrow=nv,byrow=TRUE)
      edges <- vi - this$x[1:nv,,drop=FALSE]
      n2 <- sqrt(apply(edges^2,1,sum))
      ssize <- max(c(max(n2),ssize))
    }
  }
  if (method=='sigmaplus'){
    v1 <- matrix(rep(this$x[1,],nv-1),nrow=nv-1,byrow=TRUE)
    edges <- this$x[2:nv,,drop=FALSE] - v1
    n2 <- sqrt(apply(edges^2,1,sum))
    ssize <- max(n2)
  }
  if (method=='sigmaminus'){
    v1 <- matrix(rep(this$x[1,],nv-1),nrow=nv-1,byrow=TRUE)
    edges <- this$x[2:nv,,drop=FALSE] - v1
    n2 <- sqrt(apply(edges^2,1,sum))
    ssize <- min(n2)
  }
  return(ssize)
}


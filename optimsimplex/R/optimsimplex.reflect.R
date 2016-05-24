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

optimsimplex.reflect <- function(this=NULL,fun=NULL,data=NULL){

  nv <- this$nbve
  n <- this$n
  r <- vertex() 
  r$n <- n
  r$nbve <- nv
  r$x <- matrix(0,nrow=nv,ncol=n)
  r$fv <- matrix(0,nrow=nv,ncol=1)
  r$x[1,1:n] <- this$x[1,1:n,drop=FALSE]
  r$fv[1] <- this$fv[1]
  twox1 <- matrix(rep(2*this$x[1,1:n],nv-1),nrow=nv-1,byrow=TRUE)
  r$x[2:nv,1:r$n] <- twox1 - this$x[2:nv,1:n,drop=FALSE]

  tmp <- optimsimplex.compsomefv(this=r,fun=fun,indices=2:nv,data=data)
    r <- tmp$this
    if(!is.null(data)) data <- tmp$data
  
  varargout <- list(r=r,data=data)

  return(varargout)

}


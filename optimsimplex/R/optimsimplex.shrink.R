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

optimsimplex.shrink <- function(this=NULL,fun=NULL,sigma=0.5,data=NULL){
  
  nv <- this$nbve
  mv1 <- matrix(rep(this$x[1,],nv-1),nrow=nv-1,byrow=TRUE)
  newx <- (1.0-sigma)*mv1 + sigma*this$x[2:nv,,drop=FALSE]
  this$x[2:nv,] <- newx[1:(nv-1),,drop=FALSE]

  tmp <- optimsimplex.compsomefv(this=this,fun=fun,indices=2:nv,data=data)
    this <- tmp$this
    if (!is.null(data)) data <- tmp$data

  varargout <- list(this=this,data=data)

  return(varargout)

}


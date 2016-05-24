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

vertex <- function(x,n,fv){
  
  newobj <- list(x=NULL,n=0,fv=NULL)
  
  class(newobj) <- 'vertex'

  if (!sum(missing(n),missing(x),missing(fv))%in%c(0,3))
    stop('simplex: you must provide x, n, and fv arguments or none of them.',call.=FALSE)
  
  # The dimension of the space
  if (!missing(n)){
    assert.classreal(var=n,varname='n')
    if (length(n)!=1)
      stop('simplex: The n argument must be a single value.',call.=FALSE)
    newobj$n <- n
  }
    
  # The coordinates of the vertex, with size 1 x n
  if (!missing(x)){
    assert.classreal(var=x,varname='x')
    if (size(x,1)!=1 | size(x,2)!=n)
      stop(sprintf('simplex: The x vector is expected to be a 1 x %d matrix, but current shape is %d x %d.',
                   n,size(x,1),size(x,2)),
        call.=FALSE)
    newobj$x <- x
  }
  # The function value, with size 1
  if(!missing(fv)){
    assert.classreal(var=fv,varname='fv')
    if (length(fv)!=1)
      stop('simplex: The fv argument must be a single value.',call.=FALSE)
    newobj$fv <- fv
  }
  
  return(newobj)
  
}

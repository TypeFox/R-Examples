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

simplex <- function(verbose,x,n,fv,nbve){
  
  newobj <- vertex()
  
  newobj$verbose <- FALSE 
  newobj$nbve <- 0
  
  class(newobj) <- c('simplex',class(newobj))
  
  # The verbose option, controlling the amount of messages
  if (!missing(verbose)){
    assert.classboolean(var=verbose,varname='verbose')
    newobj$verbose <- verbose
  }
  
  if (!sum(missing(n),missing(x),missing(fv),missing(nbve))%in%c(0,4))
    stop('simplex: you must provide x, n, fv, and nbve arguments or none of them.',
      call.=FALSE)
  
  # The dimension of the space
  if (!missing(n)){
    assert.classreal(var=n,varname='n')
    if (length(n)!=1)
      stop('simplex: The n argument must be a single value.',call.=FALSE)
    newobj$n <- n
  }
  
  # The coordinates of the vertices, with size nbve x n
  if (!missing(x)){
    assert.classreal(var=x,varname='x')
    if (size(x,1)!=nbve | size(x,2)!=n)
      stop(sprintf('simplex: The x argument is expected to be a %d x %d matrix, but current shape is %d x %d.',
          nbve,n,size(x,1),size(x,2)),
        call.=FALSE)
    newobj$x <- x
  }
  
  # fv: the function values, with size 1 x nbve
  if (!missing(fv)){
    assert.classreal(var=fv,varname='fv')
    if (size(fv,1)!=nbve | size(fv,2)!=1)
      stop(sprintf('simplex: The fv vector is expected to be a %d x 1 matrix, but current shape is %d x %d.',
          nbve,size(fv,1),size(fv,2)),
        call.=FALSE)
    newobj$fv <- fv
  }
  
  # The number of vertices
  if (!missing(nbve)){
    assert.classreal(var=nbve,varname='nbve')
    if (length(nbve)!=1)
      stop('simplex: The nbve argument must be a single value.',call.=FALSE)
    newobj$nbve <- nbve
  }
  
  return(newobj)

}

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

optimsimplex.setall <- function(this=NULL,simplex=NULL){

  nbve <- size(simplex,1)
  np1 <- size(simplex,2)
  if (np1>nbve)
    stop(sprintf(paste('optimsimplex_setall: The number of vertices (i.e. the number of rows)',
                       'is %d which is smaller than the number of columns %d (i.e. n+1).'),
                 nbve,np1),
         call.=FALSE)

  this$n <- np1 - 1
  this$nbve <- nbve
  this$fv[1:nbve,1] <- simplex[1:nbve,1]
  this$x[1:nbve,1:this$n] <- simplex[1:nbve,2:this$n+1]

  return(this)

}


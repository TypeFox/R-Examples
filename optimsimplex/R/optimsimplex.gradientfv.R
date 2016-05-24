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

optimsimplex.gradientfv <- function(this=NULL,fun=NULL,
                                    method='forward',data=NULL){

  if (this$nbve!=(this$n+1))
    stop(sprintf(paste('optimsimplex.gradientfv: The gradient can be applied only with ',
                       'a simplex made of n+1 points, but the dimension is %d ',
                       ' and the number of vertices is %d.',sep=''),
                 this$n,,this$nbve),
         call.=FALSE)
  if (!any(method==c('forward','centered')))
    stop(sprintf('optimsimplex.gradientfv: Unknown method %s',method),
         call.=FALSE)
  output <- list()
  if (method == 'forward'){
    g <- optimsimplex.gradforward(this=this)
  }
  if (method == 'centered'){
    tmp <- optimsimplex.gradcenter(this=this,fun=fun,data=data)
      g <- tmp$g
      if (!is.null(data)) data <- tmp$data
  }

  varargout <- list(g=g,data=data)

  return(varargout)
}


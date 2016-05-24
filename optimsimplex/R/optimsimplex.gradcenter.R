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

optimsimplex.gradcenter <- function(this=NULL,fun=NULL,data=NULL){

  g1 <- optimsimplex.gradforward(this)

  tmp <- optimsimplex.reflect(this=this,fun=fun,data=data)
    r <- tmp$r
    if (!is.null(data)) data <- tmp$data

  g2 <- optimsimplex.gradforward(r)
  g <- (g1 + g2)/2
  #r <- optimsimplex.destroy(r)

  varargout <- list(g=g,data=data)

  return(varargout)
}


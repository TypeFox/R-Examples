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
# This source code is a R port of the optimbase component
# originally written by Michael Baudin for Scilab.

optimbase.proj2bnds <- function(this=NULL,x=NULL){

  hasbounds <- optimbase.hasbounds(this=this)
  if (!hasbounds){
    p <- x
    return(p)
  }
  p <- x
  for (ix in 1:this$numberofvariables){
    xmin <- this$boundsmin[ix]
    xmax <- this$boundsmax[ix]
    pix <- p[ix]
    if (pix > xmax){
      this <- optimbase.log(this=this,
                            msg=sprintf('Projecting p(%d) = %e on max bound %e',ix,pix,xmax))
      p[ix] <- xmax
    }
    if (pix < xmin){
      this <- optimbase.log(this=this,
                            msg=sprintf('Projecting p(%d) = %e on min bound %e',ix,pix,xmin))
      p[ix] <- xmin
    }
  }

  varargout <- list(this=this,p=p)

  return(varargout)

}


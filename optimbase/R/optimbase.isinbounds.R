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

optimbase.isinbounds <- function(this=NULL,x=NULL){

  isfeasible <- TRUE
  hasbounds <- optimbase.hasbounds(this=this)
  if (hasbounds){
    for (ix in 1:this$numberofvariables){
      xmin <- this$boundsmin[ix]
      xmax <- this$boundsmax[ix]
      xix <- x[ix]
      if (xix < xmin){
        isfeasible <- FALSE
        this <- optimbase.log(this=this,
                              msg=sprintf('Component #%d/%d of x is lower than min bound %e',
                                          ix,this$numberofvariables,xmin))
        break
      }
      if (xix > xmax){
        isfeasible <- FALSE
        this <- optimbase.log(this=this,
                              msg=sprintf('Component #%d/%d of x is greater than max bound %e',
                                          ix,this$numberofvariables,xmax))
        break
      }
    }
  }

  varargout <- list(this=this,isfeasible=isfeasible)

  return(varargout)

}


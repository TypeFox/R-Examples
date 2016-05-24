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

optimbase.checkbounds <- function(this=NULL){

  maxl <- length(this$boundsmax)
  minl <- length(this$boundsmin)
  isok <- TRUE
  errmsg <- ''
  if (maxl>0 | minl>0){
    if (isok & this$numberofvariables!=maxl){
      errmsg <- sprintf(paste('The number of variables (%d) does not match the number of max bounds (%d) from [',
                              paste(this$boundsmax,collapse=' '),']\n',sep=''),
                        this$numberofvariables, maxl)
      isok <- FALSE
    }
    if (isok & this$numberofvariables!=minl){
      errmsg <- sprintf(paste('The number of variables (%d) does not match the number of min bounds (%d) from [',
                              paste(this$boundsmin,collapse=' '),']\n',sep=''),
                        this$numberofvariables, minl)
      isok <- FALSE
    }
    if (isok){
     for (ix in 1:this$numberofvariables){
#        ix <- ix+1
        xmin <- this$boundsmin[ix]
        xmax <- this$boundsmax[ix]
        if (xmax<xmin){
          errmsg <- sprintf(paste('The max bound %e for variable #%d is lower',
                                  'than the min bound %e.\n'),
                            xmax,ix,xmin)
          isok <- FALSE
          break
        }
      }
    }
  }

  varargout <- list(this=this,isok=isok,errmsg=errmsg)
  
  return(varargout)

}


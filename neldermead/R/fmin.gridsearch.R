# Copyright (C) 2010-2015 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt


fmin.gridsearch <- function(fun=NULL,x0=NULL,xmin=NULL,xmax=NULL,
    npts=3,alpha=10){
  
  gridfun <- function(x=NULL,index=NULL,fmsfundata=NULL,...){
    varargout <- list(f=fun(x),
                      g=c(),
                      c=c(),
                      gc=c(),
                      index=index,
                      this=list(costfargument=fmsfundata))
    return(varargout)
  }
  
  grid <- optimbase.gridsearch(fun=gridfun,x0=x0,xmin=xmin,xmax=xmax,
    npts=npts,alpha=alpha)
  
  return(grid)
  
}

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

optimbase.checkx0 <- function(this=NULL){

  this <- optimbase.log(this=this,msg='Checking initial guess...')
  tmp <- optimbase.isfeasible(this=this,x=this$x0)
    this <- tmp$this
    isfeasible <- tmp$isfeasible
  rm(tmp)

  isok <- (isfeasible==1)

  if (isok){
    this <- optimbase.log(this=this,msg='... initial guess is feasible.')
  } else {
    this <- optimbase.log(this=this,msg='... initial guess is not feasible.')
  }

  varargout <- list(this=this,isok=isok)

  return(varargout)

}


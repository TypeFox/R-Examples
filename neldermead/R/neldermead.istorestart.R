# Copyright (C) 2008-2009 - INRIA - Michael Baudin
# Copyright (C) 2009-2010 - DIGITEO - Michael Baudin
# Copyright (C) 2010-2015 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#
# This source code is a R port of the neldermead component
# originally written by Michael Baudin for Scilab :
# "Nelder-Mead User's Manual", 2010, Consortium Scilab - Digiteo,
# Michael Baudin, http://wiki.scilab.org/The_Nelder-Mead_Component

neldermead.istorestart <- function(this=NULL){

  status <- optimbase.get(this=this$optbase,key='status')
  if (status=='maxfuneval'){
    istorestart <- FALSE
    varargout <- list(this=this,istorestart=istorestart)
    return(varargout)
  }

  if (!any(this$restartdetection==c('oneill','kelley')))
    stop(sprintf('neldermead.istorestart: Unknown restart detection %s', this$restartdetection),
         call.=FALSE)

  if (this$restartdetection=='oneill'){
    tmp <- neldermead.isroneill(this=this)
      this <-tmp$this
      istorestart <- tmp$istorestart
    rm(tmp)
  }
  if (this$restartdetection=='kelley'){
    tmp <- neldermead.isrkelley(this=this)
      this <- tmp$this
      istorestart <- tmp$istorestart
    rm(tmp)
  }

  varargout <- list(this=this,istorestart=istorestart)

  return(varargout)

}


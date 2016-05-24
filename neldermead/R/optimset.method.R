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
# 'Nelder-Mead User's Manual', 2010, Consortium Scilab - Digiteo,
# Michael Baudin, http://wiki.scilab.org/The_Nelder-Mead_Component

optimset.method <- function(method=NULL){

  options <- list(Display=NULL,
                  FunValCheck=NULL,
                  MaxFunEvals=NULL,
                  MaxIter=NULL,
                  OutputFcn=NULL,
                  PlotFcns=NULL,
                  TolFun=NULL,
                  TolX=NULL,
                  nbMatch=NULL,
                  boundsAlpha=NULL,
                  boxScaling=NULL,
                  alphaMin=NULL)
  if (method=='fminsearch'){
    options$Display     <- 'notify'
    options$MaxFunEvals <- '200*numberofvariables'
    options$MaxIter     <- '200*numberofvariables'
    options$TolFun      <- 1e-4
    options$TolX        <- 1e-4
  } else if (method=='fminbnd'){
    options$Display     <- 'notify'
    options$MaxFunEvals <- '200*numberofvariables'
    options$MaxIter     <- '200*numberofvariables'
    options$TolFun      <- 1e-4
    options$nbMatch     <- 5
    options$boundsAlpha <- 1.e-6
    options$boxScaling  <- 0.5
    options$alphaMin    <- 1.e-5
  } else {
    stop(sprintf('optim.method: No default options available: the function \'%s\' does not exist on the path.',
                 method),
         call.=FALSE)
  }

  return(options)
  
}


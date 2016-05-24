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

neldermead.isroneill <- function(this=NULL){

  n <- optimbase.get(this=this$optbase,key='-numberofvariables')
  #
  # If required, make a vector step from the scalar step
  #
  defaultstep <- this$restartstep
  stepn <- length(defaultstep)
  if (stepn!=n){
    step <- defaultstep * transpose(rep(1,n))       # Not sure about this assignment
  } else {
    step <- defaultstep
  }
  restarteps <- this$restarteps

  x <- optimbase.get(this=this$optbase,key='xopt')
  fopt <- optimbase.get(this=this$optbase,key='fopt')
  verbose <- optimbase.get(this=this$optbase,key='verbose')

  istorestart <- FALSE
  for (ix in 1:n){
    stepix <- step[ix]
    del <- stepix * restarteps
    if (del==0.0)
      del <- .Machine$double.eps

    xix <- x[ix]
    x[ix] <- xix + del
    tmp <- optimbase.function(this=this$optbase,x=x,index=2)
      this$optbase <- tmp$this
      fv <- tmp$f
      index <- tmp$index
    rm(tmp)

    if (fv<fopt){
      istorestart <- TRUE
      if (verbose==TRUE)
        this <- neldermead.log(this=this,
                               msg=sprintf('Must restart because fv=%e at [%s] is lower than fopt=%e',
                                           fv,strvec(x),fopt))
      break
    }
    x[ix] <- xix - del
    tmp <- optimbase.function(this=this$optbase,x=x,index=2)
      this$optbase <- tmp$this
      fv <- tmp$f
      index <- tmp$index
    rm(tmp)

    if (fv<fopt){
      istorestart <- TRUE
      if (verbose==TRUE)
        this <- neldermead.log(this=this,
                               msg=sprintf('Must restart because fv=%e at [%s] is lower than fopt=%e',
                                            fv,strvec(x),fopt))
      break
    }
    x[ix] <- xix
  }

  varargout <- list(this=this,istorestart=istorestart)

  return(varargout)

}


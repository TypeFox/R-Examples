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

boxlinesearch <- function(this=NULL,n=NULL,xbar=NULL,xhigh=NULL,fhigh=NULL,rho=NULL){

  verbose <- optimbase.get(this=this$optbase,key='verbose')

  if (verbose==TRUE){
    this <- neldermead.log(this=this,msg='boxlinesearch')
    this <- neldermead.log(this=this,msg=sprintf('> xhigh=[%s], fhigh=%e',strvec(xhigh),fhigh))
    this <- neldermead.log(this=this,msg=sprintf('> xbar=[%s]',strvec(xbar)))
  }
  xr <- neldermead.interpolate(x1=xbar,x2=xhigh,fac=rho)
  if (verbose==TRUE)
    this <- neldermead.log(this=this,msg=sprintf('> xr = [%s]',strvec(xr)))
  status <- FALSE
  alphamin <- this$guinalphamin
  hasnlcons <- optimbase.hasnlcons(this=this$optbase)

  #
  # Scale from xr toward xbar until fr < fhigh and update xr
  #
  xr0 <- xr
  alpha <- 1.0
  while (alpha>alphamin){
    tmp <- optimbase.function(this=this$optbase,x=xr,index=2)
      this$optbase <- tmp$this
      fr <- tmp$f
      index <- tmp$index
      if (hasnlcons) cr <- tmp$c
    rm(tmp)

    if (fr<fhigh){
      if (verbose==TRUE)
        this <- neldermead.log(this=this,
                               msg=sprintf('fr = %e improves %e : no need for scaling for f',fr,fhigh))
      status <-TRUE
      break
    }
    alpha <- alpha * this$boxineqscaling
    if (verbose==TRUE)
      this <- neldermead.log(this=this,
                             msg=sprintf('Scaling for f with alpha=%e',alpha))
    xr <- (1.0-alpha)*xbar+alpha*xr0
    if (verbose==TRUE)
      this <- neldermead.log(this=this,msg=sprintf('> xr = %s',strvec(xr)))
  }

  # If the scaling for function improvement has failed,
  # we return.
  if (!status){
    varargout <- list(this=this,status=status,xr=xr,fr=fr)
    return(varargout)
  }

  # scaledc is set to %t if xr is updated during scaling into constraints
  # That implies that the function value is to update.
  scaledc <- FALSE

  #
  # Project xr into bounds, with an additionnal alpha inside the bounds.
  # This algo is always succesful.
  # Note:
  #   If the alpha coefficient was not used, the
  #   projectinbounds method could be used directly.
  #
  hasbounds <- optimbase.hasbounds(this=this$optbase)
  if (hasbounds){
    boxboundsalpha <- this$boxboundsalpha
    boundsmax <- optimbase.get(this=this$optbase,key='boundsmax')
    boundsmin <- optimbase.get(this=this$optbase,key='boundsmin')
    for (ix in 1:n){
      xmin <- boundsmin[ix]
      xmax <- boundsmax[ix]
      xrix <- xr[ix]
      if (xrix>xmax){
        if (verbose==TRUE)
          this <- neldermead.log(this=this,
                                 msg=sprintf('Projecting index #%d = %e on max bound %e - %e',
                                             ix,xrix,xmax,boxboundsalpha))
        xr[ix] <- xmax - boxboundsalpha
        if (!scaledc)
          scaledc <- TRUE
      }
      if (xrix<xmin){
        if (verbose==TRUE)
          this <- neldermead.log(this=this,
                                 msg=sprintf('Projecting index #%e = %e on min bound %e - %e',
                                             ix,xrix,xmin,boxboundsalpha))
        xr[ix] <- xmin + boxboundsalpha
        if (!scaledc)
          scaledc = TRUE
      }
    }
    if (verbose==TRUE)
      this <- neldermead.log(this=this,
                             msg=sprintf(' > After projection into bounds xr = [%s]',strvec(xr)))
  }

  #
  # Scale from xr to xbar into nonlinear inequality constraints
  # and update xr.
  # Set status to 0 if the process fails.
  #
  nbnlc <- optimbase.get(this=this$optbase,key='nbineqconst')
  if (nbnlc==0){
    status <-TRUE
  } else {
    status <-FALSE
    alpha <- 1.0
    xr0 <- xr
    while (alpha>alphamin){
      tmp <- optimbase.isinnonlincons(this=this$optbase,x=xr)
        this$optbase <- tmp$this
        feasible <- tmp$isfeasible
      rm(tmp)
      if (feasible){
        status <-TRUE
        break
      }
      alpha <- alpha * this$boxineqscaling
      if (verbose==TRUE)
        this <- neldermead.log(this=this,
                               msg=sprintf('Scaling for nonlinear/linear inequality constraints with alpha=%e from xbar=[%s] toward [%s]',
                                           alpha,strvec(xbar),strvec(xr0)))
      xr <- (1.0-alpha)*xbar+alpha*xr0
      if (verbose==TRUE)
        this <- neldermead.log(this=this,msg=sprintf('> xr = [%s]',strvec(xr)))
      if (!scaledc)
        scaledc <- TRUE
    }
  }

  # If scaling failed, returns immediately
  # (there is no need to update the function value).
  if (!status){
    varargout <- list(this=this,status=status,xr=xr,fr=fr)
    return(varargout)
  }
  if (scaledc){
    # Re-compute the function value at scaled point
    tmp <- optimbase.function(this=this$optbase,x=xr,index=2)
      this$optbase <- tmp$this
      fr <- tmp$f
      index <- tmp$index
      if (hasnlcons) cr <- tmp$c
    rm(tmp)
  }

  varargout <- list(this=this,status=status,xr=xr,fr=fr)

  return(varargout)

}


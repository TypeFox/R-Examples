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

scaleinconstraints <- function(this=NULL,x=NULL,xref=NULL){

  p <- x
  hasbounds <- optimbase.hasbounds(this=this$optbase)
  nbnlc <- optimbase.get(this=this$optbase,key='nbineqconst')
  varargout <- list(this=this,isscaled=NULL,p=p)
  #
  # 1. No bounds, no nonlinear inequality constraints
  # => no problem
  #
  if ((hasbounds==FALSE) & (nbnlc==0)){
    varargout$isscaled <- TRUE
    return(varargout)
  } 
  #
  # 2. Scale into bounds
  #
  if (hasbounds){
    tmp <- optimbase.proj2bnds(this=this$optbase,x=p)
      this$optbase <- tmp$this
      p <- tmp$p
    rm(tmp)
    this <- neldermead.log(this=this,
                           msg=sprintf(' > After projection into bounds p = [%s]',strvec(p)))
  }
  #
  # 3. Scale into non linear constraints
  # Try the current point and see if the constraints are satisfied.
  # If not, move the point 'halfway' to the centroid,
  # which should satisfy the constraints, if
  # the constraints are convex.
  # Perform this loop until the constraints are satisfied.
  # If all loops have been performed without success, the scaling
  # has failed.
  #
  isscaled <- FALSE
  alpha <- 1.0
  p0 <- p
  while (alpha > this$guinalphamin){
    tmp <- optimbase.isinnonlincons(this=this$optbase,x=p)
      this$optbase <- tmp$this
      feasible <- tmp$isfeasible
    rm(tmp)
    if (feasible){
      isscaled = TRUE
      break
    }
    alpha <- alpha * this$boxineqscaling
    this <- neldermead.log(this=this,
                           msg=sprintf('Scaling inequality constraint with alpha = %e',alpha))
   p <-(1.0 - alpha ) * xref + alpha * p0
  }
  this <- neldermead.log(this=this,
                         msg=sprintf(' > After scaling into inequality constraints p = [%s]',strvec(p)))
  if (!isscaled){
    this <- neldermead.log(this=this,
                           msg=sprintf(' > Impossible to scale into constraints after %d loops',
                                       this$optbase.nbineqconst))
  }

  varargout <- list(this=this, isscaled=isscaled, p=p)

  return(varargout)
  
}


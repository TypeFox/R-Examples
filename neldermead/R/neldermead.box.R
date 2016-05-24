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

neldermead.box <- function(this=this){
  
  # Check settings correspond to algo
  hascons <- optimbase.hasconstraints(this=this$optbase)
  if (!hascons)
    stop('neldermead.box: Problem has no constraints, but Box algorithm is designed for them.',
         call.=FALSE)
  verbose <- optimbase.get(this=this$optbase,key='verbose')

  #
  # Order the vertices for the first time
  #
  simplex <- this$simplex0
  n <- optimbase.get(this=this$optbase,key='numberofvariables')
  fvinitial <- optimbase.get(this=this$optbase,key='fx0')

  # Sort function values and x points by increasing function value order
  this <- neldermead.log(this=this,msg='Step #1 : order')
  simplex <- optimsimplex.sort(this=simplex)

  # Transpose, because optimsimplex returns row vectors
  currentcenter <- transpose(optimsimplex.center(this=simplex))
  currentxopt <- optimbase.get(this=this$optbase,key='x0')
  newfvmean <- optimsimplex.fvmean(this=simplex)
  nbve <- optimsimplex.getnbve(this=simplex)
  ihigh <- nbve
  inext <- ihigh - 1
  ilow <- 1
  hasbounds <- optimbase.hasbounds(this=this$optbase)
  nbnlc <- optimbase.get(this=this$optbase,key='nbineqconst')
  rho <- this$boxreflect

  #
  # Initialize
  #
  terminate <- FALSE
  iter <- 0
  step <- 'init'

  #
  # Nelder-Mead Loop
  #
  while(!terminate ){
    this$optbase <- optimbase.incriter(this=this$optbase)
    iter <- iter + 1
    xlow <- transpose(optimsimplex.getx(this=simplex,ive=ilow))
    flow <- optimsimplex.getfv(this=simplex,ive=ilow)
    xhigh <- transpose(optimsimplex.getx(this=simplex,ive=ihigh))
    fhigh <- optimsimplex.getfv(this=simplex,ive=ihigh)
    xn <- transpose(optimsimplex.getx(this=simplex,ive=inext))
    fn <- optimsimplex.getfv(this=simplex,ive=inext)
    #
    # Store history
    #
    xcoords <- optimsimplex.getallx(this=simplex)
    allfv <- optimsimplex.getallfv(this=simplex)
    this <- neldermead.storehistory(this=this,n=n,fopt=flow,xopt=xlow,fv=allfv,xcoords=xcoords)
    currentfopt <- flow
    previousxopt <- currentxopt
    currentxopt <- xlow
    previouscenter <- currentcenter
    currentcenter <- transpose(optimsimplex.center(this=simplex))
    oldfvmean <- newfvmean
    newfvmean <- optimsimplex.fvmean(this=simplex)
    if (verbose==TRUE){
      deltafv <- abs(optimsimplex.deltafvmax(this=simplex))
      totaliter <- optimbase.get(this=this$optbase,key='iterations')
      funevals <- optimbase.get(this=this$optbase,key='funevals')
      ssize <- optimsimplex.size(this=simplex)
      this <- neldermead.log(this=this,msg=sprintf('================================================================='))
      this <- neldermead.log(this=this,msg=sprintf('Iteration #%d (total = %d)',iter,totaliter))
      this <- neldermead.log(this=this,msg=sprintf('Function Eval #%d',funevals))
      this <- neldermead.log(this=this,msg=sprintf('Xopt: [%s]',strvec(xlow)))
      this <- neldermead.log(this=this,msg=sprintf('Fopt: %e',flow))
      this <- neldermead.log(this=this,msg=sprintf('DeltaFv: %e',deltafv))
      this <- neldermead.log(this=this,msg=sprintf('Center: [%s]',strvec(currentcenter)))
      this <- neldermead.log(this=this,msg=sprintf('Size: %e',ssize))
      str <- optimsimplex.tostring(x=simplex)
      for (i in 1:nbve){
        this <- neldermead.log(this=this,msg=str[i])
      }
    }
    this$optbase <- optimbase.set(this=this$optbase,key='xopt',value=xlow)
    this$optbase <- optimbase.set(this=this$optbase,key='fopt',value=flow)
    neldermead.outputcmd(this=this,state='iter',simplex=simplex,step=step)

    #
    # Update termination flag
    #
    if (iter>1){
      tmp <- neldermead.termination(this=this,fvinitial=fvinitial,oldfvmean=oldfvmean,
                                    newfvmean=newfvmean,previousxopt=previouscenter,
                                    currentxopt=currentcenter,simplex=simplex)
        this <- tmp$this
        terminate <- tmp$terminate
        status <- tmp$status
      rm(tmp)
      if (terminate) {
        if (verbose==TRUE)
          this <- neldermead.log(this=this,msg=sprintf('Terminate with status : %s',status))
        break
      }
    }

    #
    # Compute xbar, center of better vertices
    #
    if (verbose==TRUE)
      this <- neldermead.log(this=this,msg='Reflect')
    xbar <- transpose(optimsimplex.xbar(this=simplex))
    if (verbose==TRUE)
      this <- neldermead.log(this=this,msg=sprintf('xbar=[%s]',strvec(xbar)))

    #
    # Search a point improving cost function
    # and satisfying constraints.
    #
    tmp <- boxlinesearch(this=this,n=n,xbar=xbar,xhigh=xhigh,fhigh=fhigh,rho=rho)
      this <- tmp$this
      scaled <- tmp$status
      xr <- tmp$xr
      fr <- tmp$fr
    rm(tmp)
    if (!scaled){
      status <- 'impossibleimprovement'
      break
    }
    if (verbose==TRUE){
      this <- neldermead.log(this=this,msg=sprintf('xr=[%s], f(xr)=%f',strvec(xr),fr))
      this <- neldermead.log(this=this,msg='  > Perform Reflection')
    }
    simplex <- optimsimplex.setve(this=simplex,ive=ihigh,fv=fr,x=transpose(xr))
    step <- 'boxreflection'

    #
    # Sort simplex
    #
    if (verbose==TRUE)
      this <- neldermead.log(this=this,msg='Sort')
    simplex <- optimsimplex.sort(this=simplex)
  }

  this$optbase <- optimbase.set(this=this$optbase,key='xopt',value=xlow)
  this$optbase <- optimbase.set(this=this$optbase,key='fopt',value=flow)
  this$optbase <- optimbase.set(this=this$optbase,key='status',value=status)
  this$simplexopt <- simplex

  return(this)

}


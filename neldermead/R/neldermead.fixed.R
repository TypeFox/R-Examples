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

neldermead.fixed <- function(this=NULL){

  # Check settings correspond to algo
  hascons <- optimbase.hasnlcons(this$optbase)
  if (hascons)
    stop('neldermead.fixed: Problem has constraints, but fixed algorithm ignores them.',
         call.=FALSE)
  verbose <- optimbase.get(this=this$optbase,key='verbose')
  #
  # Order the vertices for the first time
  #
  simplex <- this$simplex0
  n <- optimbase.get(this=this$optbase,key='numberofvariables')
  fvinitial <- optimbase.get(this=this$optbase,key='fx0')
  # Sort function values and x points by increasing function value order
  this <- neldermead.log(this=this,msg='Sort')
  simplex <- optimsimplex.sort(this=simplex)
  #
  # Compute center of simplex
  #
  # Transpose, because optimsimplex returns row vectors
  currentcenter <- transpose(optimsimplex.center(this=simplex))
  newfvmean <- optimsimplex.fvmean(this=simplex)
  currentxopt <- optimbase.get(this=this$optbase,key='x0')
  #
  # Set indices for 'clarity'
  #
  ilow <- 1
  ihigh <- n + 1
  inext <- n
  #
  # Initialize
  #
  terminate <- FALSE
  iter <- 0
  step <- 'init'
  #
  # main N-M loop
  #
  while (!terminate){
    this$optbase <- optimbase.incriter(this=this$optbase)
    iter <- iter + 1
    xlow <- transpose(optimsimplex.getx(this=simplex,ive=ilow))
    flow <- optimsimplex.getfv(this=simplex,ive=ilow)
    xhigh <- transpose(optimsimplex.getx(this=simplex,ive=ihigh))
    fhigh <- optimsimplex.getfv(this=simplex,ive=ihigh)
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
      this <- neldermead.log(this=this,msg=sprintf('Xopt: %s',strvec(xlow)))
      this <- neldermead.log(this=this,msg=sprintf('Fopt: %e',flow))
      this <- neldermead.log(this=this,msg=sprintf('DeltaFv: %e',deltafv))
      this <- neldermead.log(this=this,msg=sprintf('Center: %s',strvec(currentcenter)))
      this <- neldermead.log(this=this,msg=sprintf('Size: %e',ssize))
      str <- optimsimplex.tostring(x=simplex)
      for (i in 1:(n+1)){
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

      if (terminate){
        if (verbose==TRUE)
          this <- neldermead.log(this=this,msg=sprintf('Terminate with status: %s',status))
        break
      }
    }
    #
    # Compute xbar, center of better vertices
    #
    if (verbose==TRUE)
      this <- neldermead.log(this=this,msg=sprintf('Reflect'))
    xbar <- transpose(optimsimplex.xbar(this=simplex))
    if (verbose==TRUE)
      this <- neldermead.log(this=this,msg=paste('xbar=',strvec(xbar),sep=''))

    #
    # Reflect the worst point with respect to center
    #
    xr <- neldermead.interpolate(x1=xbar,x2=xhigh,fac=this$rho)
    tmp <- optimbase.function(this=this$optbase,x=xr,index=2)
      this$optbase <- tmp$this
      fr <- tmp$f
      index <- tmp$index
    rm(tmp)
    if (verbose==TRUE)
      this <- neldermead.log(this=this,msg=sprintf(paste('xr=',strvec(xr),', f(xr)=%f',sep=''),fr))

    #
    # Replace worst point by xr if it is better
    #
    if (fr<fhigh){
      if (verbose==TRUE)
        this <- neldermead.log(this=this,msg='  > Perform reflect')
      simplex <- optimsimplex.setve(this=simplex,ive=ihigh,fv=fr,x=transpose(xr))
      step <- 'reflection'
    } else {
      # Reflect / xnext
      xnext <- transpose(optimsimplex.getx(this=simplex,ive=inext))
      fnext <- optimsimplex.getfv(this=simplex,ive=inext)
      xbar2 <- transpose(optimsimplex.xbar(this=simplex,iexcl=inext))
      if (verbose==TRUE)
        this <- neldermead.log(this=this,msg=paste('xbar2=',strvec(xbar2),sep=''))

      xr2 <- neldermead.interpolate(x1=xbar2,x2=xnext,fac=this$rho)
      tmp <- optimbase.function(this=this$optbase,x=xr2,index=2)
        this$optbase <- tmp$this
        fr2 <-tmp$f
        index <- tmp$index
        rm(tmp)
      if (verbose==TRUE)
        this <- neldermead.log(this=this,msg=sprintf(paste('xr2=',strvec(xr2),', f(xr2)=%f',sep=''),fr2))

      if (fr2<fnext){
        if (verbose==TRUE)
          this <- neldermead.log(this=this,msg='  > Perform reflect / next')

        simplex <- optimsimplex.setve(this=simplex,ive=inext,fv=fr2,x=transpose(xr2))
        step <- 'reflectionnext'
      } else {
        #  Shrink
        if (verbose==TRUE)
          this <- neldermead.log(this=this,msg='  > Perform Shrink')

        tmp <- optimsimplex.shrink(this=simplex,fun=costf.transposex,sigma=this$sigma,data=this)
          simplex <- tmp$this
          this <- tmp$data
        rm(tmp)
        step <- 'shrink'
      }
    }

    #
    # Sort simplex
    #
    simplex <- optimsimplex.sort(this=simplex)
  }
  this$optbase <- optimbase.set(this=this$optbase,key='xopt',value=xlow)
  this$optbase <- optimbase.set(this=this$optbase,key='fopt',value=flow)
  this$optbase <- optimbase.set(this=this$optbase,key='status',value=status)
  this$simplexopt <- simplex

  return(this)

}


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

neldermead.variable <- function(this=NULL){
  
  #
  # Check settings correspond to algo
  #
  hascons <- optimbase.hasconstraints(this=this$optbase)
  if (hascons){
    stop('neldermead.variable: Problem has constraints, but variable algorithm ignores them.',
         call.=FALSE)
  }
  verbose <- optimbase.get(this=this$optbase,key='verbose')

  #
  # Order the vertices for the first time
  #
  simplex <- this$simplex0
  n <- optimbase.get(this=this$optbase,key='numberofvariables')
  if (n==0)
    stop('neldermead.variable: The number of variable is zero.',
         call.=FALSE)
  fvinitial <- optimbase.get(this=this$optbase,key='fx0')

  # Sort function values and x points by increasing function value order
  this <- neldermead.log(this=this,msg='Step #1 : order')
  simplex <- optimsimplex.sort(this=simplex)

  # Transpose, because optimsimplex returns row vectors
  currentcenter <- transpose(optimsimplex.center(this=simplex))
  currentxopt <- optimbase.get(this=this$optbase,key='x0')
  newfvmean <- optimsimplex.fvmean(this=simplex)
  greedy <- this$greedy

  #
  # Initialize
  #
  terminate <- FALSE
  iter <- 0
  step <- 'init'

  #
  # Nelder-Mead Loop
  #
  while (!terminate){
    this$optbase <- optimbase.incriter(this=this$optbase)
    iter <- iter + 1

    # Transpose, because optimsimplex returns row vectors
    xlow <- transpose(optimsimplex.getx(this=simplex,ive=1))
    flow <- optimsimplex.getfv(this=simplex,ive=1)
    xhigh <- transpose(optimsimplex.getx(this=simplex,ive=n+1))
    fhigh <- optimsimplex.getfv(this=simplex,ive=n+1)
    xn <- transpose(optimsimplex.getx(this=simplex,ive=n))
    fn <- optimsimplex.getfv(this=simplex,ive=n)

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
        this <- neldermead.log(this=this,msg=sprintf('Terminate with status : %s',status))
        break
      }
    }

    #
    # Compute xbar, center of better vertices
    #
    if (verbose==TRUE)
      this <- neldermead.log(this=this,msg=sprintf('Reflect'))

    # Transpose, because optimsimplex returns row vectors
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
      this <- neldermead.log(this=this,
                             msg=sprintf(paste('xr=[',strvec(xr),'], f(xr)=%f',sep=''),fr))

    if (fr>=flow & fr<fn){
      # Reflexion
      if (verbose==TRUE)
        this <- neldermead.log(this=this,msg='  > Perform reflection')
      simplex <- optimsimplex.setve(this=simplex,ive=n+1,fv=fr,x=transpose(xr))
      step <- 'reflection'
    }
    if (fr<flow){
      # Expand
      if (verbose==TRUE)
        this <- neldermead.log(this=this,msg='Expand')
      xe <- neldermead.interpolate(x1=xbar,x2=xhigh,fac=this$rho*this$chi)
      tmp <- optimbase.function(this=this$optbase,x=xe,index=2)
        this$optbase <- tmp$this
        fe <- tmp$f
        index <- tmp$index
      rm(tmp)

      if (verbose==TRUE)
        this <- neldermead.log(this=this,
                               msg=sprintf(paste('xe=',strvec(xe),', f(xe)=%f',sep=''),fe))

      if (greedy){
        if (fe<flow){
          if (verbose==TRUE)
            this <- neldermead.log(this=this,msg='  > Perform Greedy Expansion')
          simplex <- optimsimplex.setve(this=simplex,ive=n+1,fv=fe,x=transpose(xe))
          step <- 'expansion'
        } else {
          if (verbose==TRUE)
            this <- neldermead.log(this=this,msg='  > Perform Greedy Reflection')
          simplex <- optimsimplex.setve(this=simplex,ive=n+1,fv=fr,x=transpose(xr))
          step <- 'reflection'
        }
      } else {
        if (fe<fr){
          if (verbose==TRUE)
            this <- neldermead.log(this=this,msg='  > Perform Expansion')
          simplex <- optimsimplex.setve(this=simplex,ive=n+1,fv=fe,x=transpose(xe))
          step <- 'expansion'
        }else{
          if (verbose==TRUE)
            this <- neldermead.log(this=this,msg='  > Perform Reflection')
          simplex <- optimsimplex.setve(this=simplex,ive=n+1,fv=fr,x=transpose(xr))
          step <- 'reflection'
        }
      }
    }
    if (fr>=fn & fr<fhigh){
      # Outside contraction
      if (verbose==TRUE)
        this <- neldermead.log(this=this,msg='Contract - outside')
      xc <- neldermead.interpolate(x1=xbar,x2=xhigh,fac=this$rho*this$gamma)
      tmp <- optimbase.function(this=this$optbase,x=xc,index=2)
        this$optbase <- tmp$this
        fc <- tmp$f
        index <- tmp$index
      rm(tmp)
      if (verbose==TRUE)
        this <- neldermead.log(this=this,
                               msg=sprintf(paste('xc=',strvec(xc),', f(xc)=%f',sep=''),fc))
      if (fc<=fr){
        if (verbose==TRUE)
          this <- neldermead.log(this=this,msg='  > Perform Outside Contraction')
        simplex <- optimsimplex.setve(this=simplex,ive=n+1,fv=fc,x=transpose(xc))
        step <- 'outsidecontraction'
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
    if (fr>=fn & fr>=fhigh){
      # Inside contraction
      if (verbose==TRUE)
        this <- neldermead.log(this=this,msg='Contract - inside')
      xc <- neldermead.interpolate(x1=xbar,x2=xhigh,fac=-this$gamma)
      tmp <- optimbase.function(this=this$optbase,x=xc,index=2)
        this$optbase <- tmp$this
        fc <- tmp$f
        index <- tmp$index
      rm(tmp)
      if (verbose==TRUE)
        this <- neldermead.log(this=this,msg=sprintf(paste('xc=',strvec(xc),', f(xc)=%f',sep=''),fc))
      if(fc < fhigh ){
        if (verbose==TRUE)
          this <- neldermead.log(this=this,msg='  > Perform Inside Contraction')
        simplex <- optimsimplex.setve(this=simplex,ive=n+1,fv=fc,transpose(xc))
        step <- 'insidecontraction'
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


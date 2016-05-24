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

fminsearch <- function(fun=NULL,x0=NULL,options=NULL,verbose=FALSE) {
  
  # Check inputs
  if (missing(fun) | is.null(fun))
    stop('fminsearch: fun argument cannot be missing or NULL.',
      call. = FALSE)
  
  if (missing(x0) | is.null(x0))
    stop('fminsearch: x0 argument cannot be missing or NULL.',
      call. = FALSE)
  
  # Set options
  defaultoptions <- optimset('fminsearch')
  if (is.null(options))
    options <- defaultoptions
  
  # Set verbose
  if (!is.logical(verbose))
    verbose <- FALSE
  verbose <- verbose[1]
  
  # Coerce x0 into a column vector
  x0 <- cbind(x0)
  
  # Compute options from the options list
  numberofvariables <- prod(size(x0))
  MaxFunEvals <- optimget(options=options,key='MaxFunEvals',value=defaultoptions$MaxFunEvals)
  MaxIter     <- optimget(options=options,key='MaxIter',value=defaultoptions$MaxIter)
  TolFun      <- optimget(options=options,key='TolFun',value=defaultoptions$TolFun)
  TolX        <- optimget(options=options,key='TolX',value=defaultoptions$TolX)
  Display     <- optimget(options=options,key='Display',value=defaultoptions$Display)
  OutputFcn   <- optimget(options=options,key='OutputFcn',value=defaultoptions$OutputFcn)
  PlotFcns    <- optimget(options=options,key='PlotFcns',value=defaultoptions$PlotFcns)
  
  # If the MaxIter option is a string, we make the assumption that it is the default 200 value.
  # If not, this is the actual value.
  if (is.character(MaxIter)){
    if (MaxIter=='200*numberofvariables'){
      MaxIter = 200 * numberofvariables
    }else{
      stop(sprintf('fminsearch: Unexpected maximum number of iterations %s.',MaxIter),
           call.=FALSE)
    }
  }
  
  # If the MaxFunEvals option is a string, this is the default 200 value
  # If not, this is the actual value.
  if (is.character(MaxFunEvals)){
    if (MaxFunEvals=='200*numberofvariables'){
      MaxFunEvals = 200 * numberofvariables
    }else{
      stop(sprintf('fminsearch: Unexpected maximum number of function evaluations %s.',
                   MaxFunEvals),
           call.=FALSE)
    }
  }
  
  # Display to shell
  if (Display=='iter'){
    cat(sprintf('%10s   %10s   %10s %17s\n','Iteration','Func-count','min f(x)','Procedure'))
  }
  
  # Prepare the data structure to pass to the output function
  fmsdata <- structure(
    list(Display=Display,
      OutputFcn=OutputFcn,
      PlotFcns=PlotFcns),
    class='optimbase.outputargs')
  
  # Prepare the data structure to pass to the cost function
  fmsfundata <- structure(
    list(Fun=fun),
    class='optimbase.functionargs')
  
  # Perform Optimization
  nm <- neldermead()
  nm <- neldermead.set(this=nm,key='x0',value=x0)
  nm <- neldermead.set(this=nm,key='numberofvariables',value=numberofvariables)
  nm <- neldermead.set(this=nm,key='simplex0method',value='pfeffer')
  nm <- neldermead.set(this=nm,key='simplex0deltausual',value=0.05)
  nm <- neldermead.set(this=nm,key='simplex0deltazero',value=0.0075)
  nm <- neldermead.set(this=nm,key='method',value='variable')
  nm <- neldermead.set(this=nm,key='function',value=fminsearch.function)
  nm <- neldermead.set(this=nm,key='costfargument',value=fmsfundata)
  nm <- neldermead.set(this=nm,key='maxiter',value=MaxIter)
  nm <- neldermead.set(this=nm,key='maxfunevals',value=MaxFunEvals)
  nm <- neldermead.set(this=nm,key='tolxmethod',value=FALSE)
  nm <- neldermead.set(this=nm,key='tolfunmethod',value=FALSE)
  nm <- neldermead.set(this=nm,key='tolssizedeltafvmethod',value=TRUE)
  nm <- neldermead.set(this=nm,key='tolsimplexizemethod',value=FALSE)
  nm <- neldermead.set(this=nm,key='toldeltafv',value=TolFun)
  nm <- neldermead.set(this=nm,key='tolsimplexizeabsolute',value=TolX)
  nm <- neldermead.set(this=nm,key='checkcostfunction',value=FALSE)
  nm <- neldermead.set(this=nm,key='outputcommand',value=fminsearch.outputfun)
  nm <- neldermead.set(this=nm,key='outputcommandarg',value=fmsdata)
  nm <- neldermead.set(this=nm,key='verbose',value=verbose)
  #nm <- neldermead.set(this=nm,key='verbosetermination',value=TRUE)
  nm <- neldermead.search(this=nm)
  fval <- neldermead.get(this=nm,key='fopt')
  status <- neldermead.get(this=nm,key='status')
  
  if (!any(status==c('maxiter','maxfuneval','tolsizedeltafv')))
    stop(sprintf('fminsearch: Unknown status %s',status),
         call.=FALSE)
  
  if (status=='maxiter'){
    if ((Display=='notify') | (Display=='iter') | (Display=='final')){
      cat(paste('fminsearch:  Exiting: Maximum number of iterations has been exceeded\n',
                '         - increase MaxIter option.\n',
                '         Current function value: ',fval,'\n',sep=''))
    }
    nm$exitflag <- FALSE
  }
  if (status=='maxfuneval'){
    if ((Display=='notify') | (Display=='iter') | (Display=='final')){
      cat(paste('fminsearch:  Exiting: Maximum number of function evaluations has been exceeded\n',
                '         - increase MaxFunEvals option.\n',
                '         Current function value: ',fval,'\n',sep=''))
    }
    nm$exitflag <- FALSE
  }
  if (status=='tolsizedeltafv'){
    nm$exitflag <- TRUE
  }
  
  nm$output <- list(algorithm ='Nelder-Mead simplex direct search',
                    funcCount =neldermead.get(this=nm,key='funevals'),
                    iterations=neldermead.get(this=nm,key='iterations'),
                    message   =sprintf('%s\n%s %e\n%s %e\n','Optimization terminated:',
                                       ' the current x satisfies the termination criteria using OPTIONS.TolX of',
                                       TolX,' and F(X) satisfies the convergence criteria using OPTIONS.TolFun of',
                                       TolFun))
  
  if ((Display=='final') | (Display=='iter')){
    if (nm$exitflag==TRUE){
      cat(nm$output$message)
    }
  }
  
  return(nm)
  
}

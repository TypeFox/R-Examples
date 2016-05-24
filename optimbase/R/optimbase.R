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

optimbase <- function(verbose, x0, fx0, xopt, fopt, tolfunabsolute, 
  tolfunrelative, tolfunmethod, tolxabsolute, tolxrelative, tolxmethod, 
  maxfunevals, funevals, maxiter, iterations, fun, status, historyxopt,
  historyfopt, verbosetermination, outputcommand, outputcommandarg,
  numberofvariables, storehistory, costfargument, boundsmin, boundsmax,
  nbineqconst, logfile, logfilehandle, logstartup, withderivatives){
  
  newobj <- list(
    verbose = FALSE,
    x0 = c(),
    fx0 = c(),
    xopt = 0,
    fopt = 0,
    tolfunabsolute = 0,
    tolfunrelative = .Machine$double.eps,
    tolfunmethod = FALSE,
    tolxabsolute = 0,
    tolxrelative = .Machine$double.eps,
    tolxmethod = TRUE,
    maxfunevals = 100,
    funevals = 0,
    maxiter = 100,
    iterations = 0,
    fun = '',
    status = '',
    historyxopt = c(),
    historyfopt = c(),
    verbosetermination = FALSE,
    outputcommand = '',
    outputcommandarg = '',
    numberofvariables = 0,
    storehistory = FALSE,
    costfargument = '',
    boundsmin = c(),
    boundsmax = c(),
    nbineqconst = 0,
    logfile = '',
    logfilehandle = 0,
    logstartup = FALSE,
    withderivatives = FALSE)
  
  class(newobj) <- 'optimbase'
  
  # The verbose option, controlling the amount of messages
  if (!missing(verbose))
    newobj <- optimbase.set(this=newobj,key='verbose',value=verbose)
  
  # The initial guess
  if (!missing(x0))
    newobj <- optimbase.set(this=newobj,key='x0',value=x0)
  
  # The value of the function for the initial guess
  if (!missing(fx0))
    newobj <- optimbase.set(this=newobj,key='fx0',value=fx0)
  
  # The optimum parameter
  if (!missing(xopt))
    newobj <- optimbase.set(this=newobj,key='xopt',value=xopt)
  
  # The optimum function value
  if (!missing(fopt))
    newobj <- optimbase.set(this=newobj,key='fopt',value=fopt)
  
  # The absolute tolerance on function value
  if (!missing(tolfunabsolute))
    newobj <- optimbase.set(this=newobj,key='tolfunabsolute',value=tolfunabsolute)
  
  # The relative tolerance on function value
  if (!missing(tolfunrelative))
    newobj <- optimbase.set(this=newobj,key='tolfunrelative',value=tolfunrelative)
  
  # Possible values : FALSE, TRUE
  # This criteria is suitable for functions which minimum is
  # associated with a function value equal to 0.
  if (!missing(tolfunmethod))
    newobj <- optimbase.set(this=newobj,key='tolfunmethod',value=tolfunmethod)
  
  # The absolute tolerance on x
  if (!missing(tolxabsolute))
    newobj <- optimbase.set(this=newobj,key='tolxabsolute',value=tolxabsolute)
  
  # The relative tolerance on x
  if (!missing(tolxrelative))
    newobj <- optimbase.set(this=newobj,key='tolxrelative',value=tolxrelative)
  
  # Possible values : FALSE, TRUE
  if (!missing(tolxmethod))
    newobj <- optimbase.set(this=newobj,key='tolxmethod',value=tolxmethod)
  
  # The maximum number of function evaluations
  if (!missing(maxfunevals))
    newobj <- optimbase.set(this=newobj,key='maxfunevals',value=maxfunevals)
  
  # The number of function evaluations
  if (!missing(funevals))
    newobj <- optimbase.set(this=newobj,key='funevals',value=funevals)
  
  # The maximum number of iterations
  if (!missing(maxiter))
    newobj <- optimbase.set(this=newobj,key='maxiter',value=maxiter)
  
  # The number of iterations
  if (!missing(iterations))
    newobj <- optimbase.set(this=newobj,key='iterations',value=iterations)
  
  # The cost function
  if (!missing(fun))
    newobj <- optimbase.set(this=newobj,key='fun',value=fun)
  
  # The status of the optimization
  if (!missing(status))
    newobj <- optimbase.set(this=newobj,key='status',value=status)
  
  # The array to store the history for fopt
  if (!missing(historyfopt))
    newobj <- optimbase.set(this=newobj,key='historyfopt',value=historyfopt)
  
  # The array to store the history for xopt
  if (!missing(historyxopt))
    newobj <- optimbase.set(this=newobj,key='historyxopt',value=historyxopt)
  
  # The verbose option for termination criteria
  if (!missing(verbosetermination))
    newobj <- optimbase.set(this=newobj,key='verbosetermination',value=verbosetermination)
  
  # The command called back for output
  if (!missing(outputcommand))
    newobj <- optimbase.set(this=newobj,key='outputcommand',value=outputcommand)
  
  # The outputcommand argument is initialized as a string.
  # If the user configures this option, it is expected
  # that a matrix of values or a list, tlist, mlist is
  # passed so that the argument is appended to the name of the
  # function.
  if (!missing(outputcommandarg))
    newobj <- optimbase.set(this=newobj,key='outputcommandarg',value=outputcommandarg)
  
  # The number of variables to optimize
  if (!missing(numberofvariables))
    newobj <- optimbase.set(this=newobj,key='numberofvariables',value=numberofvariables)
  
  # The flag which enables/disables the storing of the history
  if (!missing(storehistory))
    newobj <- optimbase.set(this=newobj,key='storehistory',value=storehistory)
  
  # The costf argument is initialized as a string.
  # If the user configures this option, it is expected
  # that a matrix of values or a list, tlist, mlist is
  # passed so that the argument is appended to the name of the
  # function.
  if (!missing(costfargument))
    newobj <- optimbase.set(this=newobj,key='costfargument',value=costfargument)
  
  # Minimum bounds for the parameters
  if (!missing(boundsmin))
    newobj <- optimbase.set(this=newobj,key='boundsmin',value=boundsmin)
  
  # Maximum bounds for the parameters
  if (!missing(boundsmax))
    newobj <- optimbase.set(this=newobj,key='boundsmax',value=boundsmax)
  
  # The number of nonlinear inequality constraints
  if (!missing(nbineqconst))
    newobj <- optimbase.set(this=newobj,key='nbineqconst',value=nbineqconst)
  
  # The name of the log file
  if (!missing(logfile))
    newobj <- optimbase.set(this=newobj,key='logfile',value=logfile)
  
  # The handle for the log file
  if (!missing(logfilehandle))
    newobj <- optimbase.set(this=newobj,key='logfilehandle',value=logfilehandle)
  
  # Set to TRUE when the logging is started up
  if (!missing(logstartup))
    newobj <- optimbase.set(this=newobj,key='logstartup',value=logstartup)
  
  # Set to TRUE when the method uses derivatives
  if (!missing(withderivatives))
    newobj <- optimbase.set(this=newobj,key='withderivatives',value=withderivatives)
  
  class(newobj) <- 'optimbase'
  
  return(newobj)
  
}

optimbase.outputargs <- function(...){
  
  structure(list(...), class='optimbase.outputargs')
  
}

optimbase.functionargs <- function(...){
  
  structure(list(...), class='optimbase.functionargs')
  
}

# File:   snomadr.R
# Author: Zhenghua Nie
# Date:   Mon 16 May 2011
#
# We use ipoptr developed by Jelmer Ypma as the prototype of this package.
# Some code is copied and edited from ipoptr.
# Please reference the license of ipoptr.
#
# Input:
#   n  : the number of variables
#    x0 : vector with initial values
#    eval.f : function to evaluate objective function
#    lb : lower bounds of the control
#    ub : upper bounds of the control
#    opts : list with options that are passed to Ipopt
#       ... : arguments that will be passed to user-defined functions
#
# Output: structure with inputs and
#    call : the call that was made to solve
#    status : integer value with the status of the optimization (0 is success)
#    message : more informative message with the status of the optimization
#    bbe : number of the objective function that were executed
#    iterations : number of iterations that were executed
#    objective : value if the objective function in the solution
#    solution : optimal value of the controls
#
# Copyright (C) 2011 Zhenghua Nie. All Rights Reserved.
# This code is published under GNU GENERAL PUBLIC LICENSE.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License,  or
# (at your option) any later version.
#
# This program is distributed WITHOUT ANY WARRANTY. See the
# GNU General Public License for more details.
#  
# If you do not have a copy of the GNU General Public License,
# write to the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

snomadr <-
function( eval.f,
          n,
          bbin = NULL,
          bbout = NULL,
          x0 = NULL,
          lb = NULL,
          ub = NULL,
          nmulti = 0,   #0: call single nomad,
          random.seed = 0,  # seed will be used for generating multiple initial points.
          opts = list(),
          print.output = TRUE,  #0: if it is FALSE,  there will be no output in snomadr, if DISPLAY_DEGREE=0 and print_output is true, there will be no any output.
          information = list(),
          snomadr.environment = new.env(),
          ... ) {

  ## Save seed prior to setting

  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }

  set.seed(random.seed)

  if(length(information) > 0) {
    sinfo <- NULL
    sversion <- NULL
    shelp <- NULL
    if(!is.null(information$info)) sinfo <-information$info
    if(!is.null(information$version)) sversion <- information$version
    if(!is.null(information$help)) shelp <- information$help

    ret <- list( info=sinfo, version=sversion, help=shelp, snomadr.environment)

    attr(ret, "class") <- "snomadr"

    ## Add the current call to the list
    ret$call <- match.call()

    ## Pass snomadr object to C code
    solution <- .Call( snomadRInfo, ret )
    ## Remove the environment from the return object
    ret$environment <- NULL
    ## We have not implemented the following output from snomadRInfo
    ##ret$Info<-solution$Info
    ##ret$Version<-solution$Version
    ##ret$Help<-solution$Help

    return(ret)

  }

  information <- NULL  ## We will check whether it is NULL in print.snomadr.

  ## The number of variables should not be null.
  if (missing(n ) || missing(eval.f)) stop("Must provide the objective function and the number of variables")
  if(missing(nmulti)||nmulti < 0) nmulti <- 0
  if(missing(print.output)) print.output <- TRUE

  ## Define 'continuous' to types of variables
  if (is.null(bbin) ) { bbin <- rep (0, n)}
  ## Define 'NOMAD::OBJ' to the output type of the function
  if (is.null(bbout)) {bbout <- rep(0, 1)}

  ## Define 'infinite' lower and upper bounds of the control if they haven't been set
  if ( is.null( lb ) ) { lb <- rep( -Inf, n ) }
  if ( is.null( ub ) ) { ub <- rep(  Inf, n ) }

  ## We don't need to generate the initial point for multiple mads runs.
  if(is.null(x0)&&nmulti < 1){
    x0<-rep(0.0, n)
    for(i in 1:n){
      x0[i] <- runif(1, min=lb[i], max=ub[i])
    }
  }

  ## Change the environment of the functions that we're calling the
  ## environment of the eval.f is changed below (if it exists)
  environment( eval.f ) <- snomadr.environment

  ## Internal function to check the arguments of the functions
  checkFunctionArguments <- function( fun, arglist, funname ) {
    if( !is.function(fun) ) stop(paste(funname, " must be a function\n", sep = ""))

    ## Determine function arguments
    fargs <- formals(fun)

    if ( length(fargs) > 1 ) {
      ## Determine argument names user-defined function
      argnames.udf <- names(fargs)[2:length(fargs)]  ## remove first argument, which is x

      ## Determine argument names that where supplied to snomadr()
      argnames.supplied <- names(arglist)

      ## Determine which arguments where required but not supplied
      m1 = match(argnames.udf, argnames.supplied)
      if( any(is.na(m1)) ){
        mx1 = which( is.na(m1) )
        for( i in 1:length(mx1) ){
          stop(paste(funname, " requires argument '", argnames.udf[mx1], "' but this has not been passed to the 'snomadr' function.\n", sep = ""))
        }
      }

      ## Determine which arguments where supplied but not required
      m2 = match(argnames.supplied, argnames.udf)
      if( any(is.na(m2)) ){
        mx2 = which( is.na(m2) )
        for( i in 1:length(mx2) ){
          stop(paste("'", argnames.supplied[mx2], "' passed to (...) in 'snomadr' but this is not required in the ", funname, " function.\n", sep = ""))
        }
      }
    }
    return( 0 )
  }

  ## Extract list of additional arguments and check user-defined
  ## functions There is an error when building the package crs, so
  ## Zhenghua commented the the following two lines.

  arglist <- list(...)
  checkFunctionArguments( eval.f, arglist, 'eval.f' )

  ## Write wrappers around user-defined functions to pass additional
  ## arguments
  eval.f.wrapper <- function(x){ eval.f(x,...) }

  ## Build snomadr object
  ret <- list("eval.f"=eval.f.wrapper,
              "n"=as.integer(n),
              "bbin"=as.integer(bbin),
              "bbout"=as.integer(bbout),
              "x0"=x0,
              "lower.bounds"=lb,
              "upper.bounds"=ub,
              "nmulti"=as.integer(nmulti),
              "random.seed"=as.integer(random.seed),
              "options"=get.option.types(opts),
              "print.output"=print.output,
              "snomadr.environment"=snomadr.environment)

  attr(ret, "class") <- "snomadr"

  ## Add the current call to the list
  ret$call <- match.call()

  ## Pass snomadr object to C code
  if(nmulti == 0){
    solution <- .Call( snomadRSolve, ret )
  } else {
    solution <- .Call( smultinomadRSolve, ret )
  }

  ## Remove the environment from the return object
  ret$environment <- NULL

  ## Add solution variables to object
  ret$status <- solution$status
  ret$message <- solution$message
  ret$bbe <- solution$bbe
  ret$iterations <- solution$iterations
  ret$objective <- solution$objective
  ret$solution <- solution$solution

  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)

  return( ret )

}

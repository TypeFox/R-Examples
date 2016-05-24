## initialize methods and constructors (generating functions)
## the initialize method is called by new
## simObj initialize first calls the default
## initialize and then (if available) an object (instance!) specific
## function stored in slot "initfunc" where "initfunc" has access to all
## functions stored in the equations list (see ?addtoenv)
##
## the constructors are provided for compatibility with older releases
## of simecol and may / or may not be deprecated in future versions

setMethod("initialize", signature(.Object="simObj"),
  function(.Object, ...) {
    .Object <- callNextMethod()
    ## delete outputs from former simulations
    .Object@out <- NULL
    if (is.function(.Object@initfunc)) {
      initfunc                <- .Object@initfunc
      equations               <- .Object@equations
      environment(initfunc)   <- environment()
      environment(main)       <- environment()
      #attach(equations)
      #on.exit(detach(equations))
      equations               <- addtoenv(equations)
      .Object                 <- initfunc(.Object)
    }
    invisible(.Object)
  }
)

## constructors ('generating functions')
simObj  <- function(obj, main=NULL, equations=NULL,
                      times=c(from=0, to=10, by=1),
                      init=matrix(0), parms=list(), inputs=NULL,
                      solver="iteration", initfunc=NULL) {
  if (is(obj, "simObj")) {
    obj <- initialize(obj)
  } else {
    obj <- new("simObj", main=main,
               equations=equations, times=times,
               init=init, parms=parms, inputs=inputs, solver=solver,
               initfunc=initfunc)
  }
  invisible(obj)
}

odeModel  <- function(obj=NULL, main=NULL, equations=NULL,
                      times=c(from=0, to=10, by=1),
                      init=numeric(0), parms=numeric(0), inputs=NULL,
                      solver="rk4", initfunc=NULL) {
                      
  if (is(obj, "odeModel")) {
    obj <- initialize(obj)
  } else {
    obj <- new("odeModel", main=main,
               equations=equations, times=times,
               init=init, parms=parms, inputs=inputs, solver=solver,
               initfunc=initfunc)
  }
  invisible(obj)
}

gridModel  <- function(obj=NULL, main=NULL, equations=NULL,
                      times=c(from=0, to=10, by=1),
                      init=matrix(0), parms=list(), inputs=NULL,
                      solver="iteration", initfunc=NULL) {
  if (is(obj, "gridModel")) {
    obj <- initialize(obj)
  } else {
    obj <- new("gridModel", main=main,
               equations=equations, times=times,
               init=init, parms=parms, inputs=inputs, solver=solver,
               initfunc=initfunc)
  }
  invisible(obj)
}

rwalkModel  <- function(obj=NULL, main=NULL, equations=NULL,
                      times=c(from=0, to=10, by=1),
                      init=NULL, parms=list(), inputs=NULL,
                      solver="iteration", initfunc=NULL) {
  if (is(obj, "rwalkModel")) {
    obj <- initialize(obj)
  } else {
    obj <- new("rwalkModel", main=main,
               equations=equations, times=times,
               init=init, parms=parms, inputs=inputs, solver=solver,
               initfunc=initfunc)
  }
  invisible(obj)
}

indbasedModel  <- function(obj=NULL, main=NULL, equations=NULL,
                      times=c(from=0, to=10, by=1),
                      init=NULL, parms=list(), inputs=NULL,
                      solver="iteration", initfunc=NULL) {
  if (is(obj, "indbasedModel")) {
    obj <- initialize(obj)
  } else {
    obj <- new("indbasedModel", main=main,
               equations=equations, times=times,
               init=init, parms=parms, inputs=inputs, solver=solver,
               initfunc=initfunc)
  }
  invisible(obj)
}



## === template to derive your own initialize method ===
#setMethod("initialize",
#          signature(.Object="odeModel"),
#          function(.Object, ...) {
#            ## ~~ put your initialization code here ~~
#            callNextMethod()
#          }
#)

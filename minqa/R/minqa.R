##' Handle the common arguments to all the minimizers
##'
##' Establish defaults for the elements of the control list and parse
##' the list that was provided.  The list is converted to an
##' environment.
##'
##' @param ctrl list of control settings
##' @param n length of the par vector
##' 
##' @return an environment containing the control settings
commonArgs <- function(par, fn, ctrl, rho) {
    rho$n <- n <- length(rho$par <- as.double(par))
    stopifnot(all(is.finite(par)),
              is.function(fn),
              length(formals(fn)) >= 1)
    rho$.feval. <- integer(1)           # function evaluation counter

    ## We use all possible control settings in the default.
    ## Extra control settings are ignored.
##    cc <- do.call(function(npt = min(n+6L, 2L * n + 1L), rhobeg = NA,
##                           rhoend = NA, iprint = 0L, maxfun=10000L,
##                           obstop=TRUE, force.start=FALSE)
    cc <- do.call(function(npt = min(n+2L, 2L * n), rhobeg = NA,
                           rhoend = NA, iprint = 0L, maxfun=10000L,
                           obstop=TRUE, force.start=FALSE,...) {

      if (length(list(...))>0) warning("unused control arguments ignored")
      list(npt = npt, rhobeg = rhobeg, rhoend = rhoend,
           iprint = iprint, maxfun = maxfun, obstop = obstop,
           force.start = force.start)
      }, ctrl)

    ## Create and populate an environment
    ctrl <- new.env(parent = emptyenv()) # ctrl environment should not chain
    lapply(names(cc), function(nm) assign(nm, cc[[nm]], envir = ctrl))

    ## Adjust and check npt
    ctrl$npt <- as.integer(max(n + 2L, min(ctrl$npt, ((n+1L)*(n+2L)) %/% 2L)))
    if (ctrl$npt > (2 * n + 1))
        warning("Setting npt > 2 * length(par) + 1 is not recommended.")

    ## Check and adjust rhobeg and rhoend
    if (is.na(ctrl$rhobeg))
        ctrl$rhobeg <- min(0.95, 0.2 * max(abs(par)))
    if (is.na(ctrl$rhoend)) ctrl$rhoend <- 1.0e-6 * ctrl$rhobeg
    stopifnot(0 < ctrl$rhoend, ctrl$rhoend <= ctrl$rhobeg)

    ## Check recommended range of maxfun
    if (ctrl$maxfun < 10 * n^2)
        warning("maxfun < 10 * length(par)^2 is not recommended.")
    ctrl
}

##' Nonlinear optimization with box constraints
##'
##' Minimize a function of many variables subject to box constraints
##' by a trust region method that forms quadratic models by 
##' interpolation, using the BOBYQA code written by Mike Powell.
##' 
##' @param par numeric vector of starting parameters (length > 1)
##' @param fn function to be minimized.  The first argument must be
##'     the parameters.
##' @param lower a numeric vector of lower bounds.  If of length 1 it
##'     is expanded.
##' @param upper a numeric vector of upper bounds.  Also may be scalar.
##' @param control a list of control settings
##' @param ... optional, additional arguments to fn
##'
##' @return a list with S3 class bobyqa
##' 

bobyqa <- function(par, fn, lower = -Inf, upper = Inf, control = list(), ...)
{
    nn <- names(par) 
    ctrl <- commonArgs(par, fn, control, environment())
    n <- length(par)
    fn1 <- function(x) {  # fn1 takes exactly 1 argument
      names(x) <- nn 
      fn(x, ...) 
    }
    checkObj <- fn1(par)
    if(length(checkObj) > 1 || !is.numeric(checkObj))
      stop("Objective function must return a single numeric value.")
    ## check the upper and lower arguments, adjusting if necessary
    lower <- as.double(lower); upper <- as.double(upper)
    if (length(lower) == 1) lower <- rep(lower, n)
    if (length(upper) == 1) upper <- rep(upper, n)
    stopifnot(length(lower) == n, length(upper) == n, all(lower < upper))
    if (any(par < lower | par > upper)) {
        if (ctrl$obstop)
            stop("Starting values violate bounds")
        else {
            par <- pmax(lower, pmax(par, upper))
            warning("Some parameters adjusted to nearest bound")
        }
    }
    rng <- upper - lower

    if (any(rng < 2 * ctrl$rhobeg)) {
        warning("All upper - lower must be >= 2*rhobeg. Changing rhobeg") 
        ctrl$rhobeg <- 0.2 * min(rng)
    }
    
    verb <- 1 < (ctrl$iprint <- as.integer(ctrl$iprint))
    ## Modifications to par if too close to boundary
    if (all(is.finite(upper)) && all(is.finite(lower)) &&
        all(par >= lower) && all(par <= upper) ) {
      if (verb) cat("ctrl$force.start = ", ctrl$force.start,"\n")
      if (!ctrl$force.start) {
        i <- rng < ctrl$rhobeg # Jens modification
        if (any(i)) {
          par[i] <- lower[i] + ctrl$rhobeg
          warning("Some parameters adjusted away from lower bound")
        }
        i <- rng < ctrl$rhobeg      # Jens modification
        if (any(i)) {
                par[i] <- upper[i] - ctrl$rhobeg
                warning("Some parameters adjusted away from upper bound")
              }
      }
    }
    if (verb) {
      cat("npt =", ctrl$npt, ", n = ",n,"\n")
      cat("rhobeg = ", ctrl$rhobeg,", rhoend = ", ctrl$rhoend, "\n")
      
    }
    if(ctrl$iprint > 0)
      cat("start par. = ", par, "fn = ", checkObj, "\n")      
        
    retlst<- .Call(bobyqa_cpp, par, lower, upper, ctrl, fn1)
# JN 20100810 
    if (retlst$ierr > 0){
##	cat("ierr = ",retlst$ierr,"\n")
##        newuoa allowed ierr in c(10, 20, 320, 390, 430)
  	if (retlst$ierr == 10) {
		retlst$ierr<-2
		retlst$msg<-"bobyqa -- NPT is not in the required interval"
	} else if (retlst$ierr == 320) {
		retlst$ierr<-5
		retlst$msg<-"bobyqa detected too much cancellation in denominator"
	} else if (retlst$ierr == 390) {
		retlst$ierr<-1
		retlst$msg<-"bobyqa -- maximum number of function evaluations exceeded"
	} else if (retlst$ierr == 430) {
		retlst$ierr<-3
		retlst$msg<-"bobyqa -- a trust region step failed to reduce q"
	} else if (retlst$ierr == 20) {
		retlst$ierr<-4
		retlst$msg<-"bobyqa -- one of the box constraint ranges is too small (< 2*RHOBEG)"
        }
    } else { 
	retlst$msg<-"Normal exit from bobyqa"
    }
    retlst # return(retlst)
}

##' An R interface to the NEWUOA implementation of Powell
##'
##' Minimize a function of many variables by a trust region method
##' that forms quadratic models by interpolation, using the NEWUOA
##' code written by Mike Powell.
##' 
##' @param par numeric vector of starting parameters (length > 1)
##' @param fn function to be minimized.  The first argument must be
##'     the parameters.
##' @param control a list of control settings
##' @param ... optional, additional arguments to fn
##'
##' @return a list with S3 class c("newuoa", "minqa")
newuoa <- function(par, fn, control = list(), ...)
{
    nn <- names(par) 
    ctrl <- commonArgs(par + 0, fn, control, environment())
    n <- length(par)
    fn1 <- function(x) {  # fn1 takes exactly 1 argument
      names(x) <- nn 
      fn(x, ...) 
    }
    checkObj <- fn1(par)
    if(length(checkObj) > 1 || !is.numeric(checkObj))
      stop("Objective function must return a single numeric value.")
    verb <- 1 < (ctrl$iprint <- as.integer(ctrl$iprint))
    if (verb) {
      cat("npt =", ctrl$npt, ", n = ",n,"\n")
      cat("rhobeg = ", ctrl$rhobeg,", rhoend = ", ctrl$rhoend, "\n")
    }
    if(ctrl$iprint > 0)
      cat("start par. = ", par, "fn = ", checkObj, "\n")
      
    retlst<-.Call(newuoa_cpp, par, ctrl, fn1)
# JN 20100810 
    if (retlst$ierr > 0){
##	cat("ierr = ",retlst$ierr,"\n")
##        newuoa allowed ierr in c(10, 320, 390, 3701)
  	if (retlst$ierr == 10) {
		retlst$ierr<-2
		retlst$msg<-"newuoa -- NPT is not in the required interval"
	} else if (retlst$ierr == 320) {
		retlst$ierr<-5
		retlst$msg<-"newuoa detected too much cancellation in denominator"
	} else if (retlst$ierr == 390) {
		retlst$ierr<-1
		retlst$msg<-"newuoa -- maximum number of function evaluations exceeded"
	} else if (retlst$ierr == 3701) {
		retlst$ierr<-3
		retlst$msg<-"newuoa -- a trust region step failed to reduce q"
        }
    } else { 
	retlst$msg<-"Normal exit from newuoa"
    }
    retlst # return(retlst)
}




##' An R interface to the UOBYQA implementation of Powell
##'
##' Minimize a function of many variables by a trust region method
##' that forms quadratic models by interpolation, using the UOBYQA
##' code written by Mike Powell.
##' 
##' @param par numeric vector of starting parameters (length > 1)
##' @param fn function to be minimized.  The first argument must be
##'     the parameters.
##' @param control a list of control settings
##' @param ... optional, additional arguments to fn
##'
##' @return a list with S3 class uobyqa
uobyqa <- function(par, fn, control = list(), ...)
{   nn <- names(par) 
    ctrl <- commonArgs(par + 0, fn, control, environment())
    n <- length(par)
    fn1 <- function(x) {  # fn1 takes exactly 1 argument
      names(x) <- nn 
      fn(x, ...) 
    }
    checkObj <- fn1(par)
    if(length(checkObj) > 1 || !is.numeric(checkObj))
      stop("Objective function must return a single numeric value.")
    verb <- 1 < (ctrl$iprint <- as.integer(ctrl$iprint))
    if (verb) {
      cat("npt =", ctrl$npt, ", n = ",n,"\n")
      cat("rhobeg = ", ctrl$rhobeg,", rhoend = ", ctrl$rhoend, "\n")
    }
    if(ctrl$iprint > 0)
      cat("start par. = ", par, "fn = ", checkObj, "\n")
      
    retlst<-.Call(uobyqa_cpp, par, ctrl, fn1)
    
# JN 20100810 
    if (retlst$ierr > 0){
##	cat("ierr = ",retlst$ierr,"\n")
##       uobyqa allowed ierr in c(390, 2101)
	if (retlst$ierr == 390) {
		retlst$ierr<-1
		retlst$msg<-"uobyqa -- maximum number of function evaluations exceeded"
	} else if (retlst$ierr == 2101) {
		retlst$ierr<-3
		retlst$msg<-"uobyqa -- a trust region step failed to reduce q"
	}
    } else { 
	retlst$msg<-"Normal exit from uobyqa"
    }
    retlst # return(retlst)
}

##' Print method for minqa objects (S3)
##'
##' @param x an object of class that inherits from minqa
##' @param digits number of significant digits - doesn't seem to be used
##' @param ... optional arguments.  None are used.
##' @return invisible(x) - side effect is to print
##' @author Douglas Bates
print.minqa <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("parameter estimates:", toString(x$par), "\n")
  cat("objective:", toString(x$fval), "\n")
  cat("number of function evaluations:", toString(x$feval), "\n")
  
  invisible(x)
}

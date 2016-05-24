### ============================================================================
###
### timelags and delay differential equations
###
### ============================================================================

## =============================================================================
## lagged values and derivates are obtained in the R-code via functions 
## lagvalue and lagderiv 
## =============================================================================

lagvalue <- function (t, nr=NULL) {
  if (is.null(nr)) nr <- 0
  out <- .Call("getLagValue", t = t, PACKAGE = "deSolve", as.integer(nr))
  return(out)
}

lagderiv <- function (t, nr=NULL) {
  if (is.null(nr)) nr <- 0
  out <- .Call("getLagDeriv", t = t, PACKAGE = "deSolve", as.integer(nr))
  return(out)
}

### ============================================================================
### solving Delay Differential Equations
### ============================================================================

dede <- function(y, times, func=NULL, parms, method = c( "lsoda", "lsode", 
    "lsodes", "lsodar", "vode", "daspk", "bdf", "adams", "impAdams", "radau"),
    control=NULL,  ...) {
    if (is.null(control)) control <- list(mxhist = 1e4)

    if (is.null(method)) 
        method <- "lsoda"
    else if (is.function(method)) 
        res <- method(y, times, func, parms, lags = control, ...)
    else if (is.complex(y))
     stop ("cannot run dede with complex y")
    else 
      res <- switch(match.arg(method), 
       lsoda = lsoda(y, times, func, parms, lags = control, ...), 
       vode = vode(y, times, func, parms, lags = control, ...), 
       lsode = lsode(y, times, func, parms, lags = control, ...), 
       lsodes = lsodes(y, times, func, parms, lags = control, ...), 
       lsodar = lsodar(y, times, func, parms, lags = control, ...), 
       daspk = daspk(y, times, func, parms, lags = control, ...),
       bdf  = lsode(y, times, func, parms, mf = 22, lags = control, ...),
       adams = lsode(y, times, func, parms, mf = 10, lags = control, ...), 
       radau = radau(y, times, func, parms, lags = control, ...),
       impAdams = lsode(y, times, func, parms, mf = 12, lags = control, ...)
    )
    return(res)
}


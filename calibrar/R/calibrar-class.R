
# Methods for calibrar.results class ---------------------------------------

#' @export
print.calibrar.results = function(x, ...) {
  
  cat("Calibration finished.\n")
  cat("Function value:", x$value, "\n")
  cat("Parameters:\n")
  print(x=x$par, ...)
  if(!all(x$active)) cat("* Parameters not calibrated.\n")
   
}

#' @export
coef.calibrar.results = function(object, ...) {
  return(object$par)
}

#' @export
predict.calibrar.results = function(object, ...) {

  obj = object$fn
  
  if(is.null(attr(obj, which="fn"))) {
    warning("Predict is only available for functions created with 'createObjectiveFunction'.")
    return(invisible())
  }
  
  fn = match.fun(attr(obj, which="fn"))  
  out = fn(object$par)
  return(out)
}

#' @export
plot.calibrar.results = function(x, ...) {
  return(invisible(NULL))
}

#' @export
summary.calibrar.results = function(object, ..., pars=NULL) {
  oNames = as.character(match.call())[-1]
  objs = list(...)
  useDots = all(sapply(objs, FUN=inherits, what="calibrar.results"))
  if(useDots & length(objs)>0) {
    objs = c(list(object), objs)
    out = lapply(objs, FUN=.summaryCalibrarResults, pars=pars)
    out = do.call(rbind, out)
    rownames(out) = oNames[seq_along(objs)]
    class(out) = c("summary.calibrar.results", class(out))
    return(out)
  }
  object$nphases = length(object$phases)
  object$nactive = sum(object$active)
  object$npar = length(unlist(object$par))
  class(object) = "summary.calibrar.results"
  return(object)
}

.summaryCalibrarResults = function(x, pars=NULL) {
  ox = x
  xpars = unclass(unlist(x$par))
  npar = length(xpars)
  if(is.null(pars)) pars = seq_len(npar)
  if(is.character(pars)) {
    if(all(pars %in% names(x$par))) {
      xpars = xpars[pars]
      xpars = unclass(unlist(xpars))
    }
  } else {
    if(all(pars %in% seq_len(npar)))
      xpars = unlist(unclass(xpars))[pars]
  }
  
  out = c(value=x$value, xpars)
  return(out)
}

#' @export
print.summary.calibrar.results = function(x, digits=3, ...) {
  if(is.matrix(x)) {
    print(unclass(x), digits=digits, ...)
    return(invisible())
  }
  cat(sprintf("Calibration in %d %s.\n", x$nphases, 
              ifelse(x$nphases==1, "phase", "phases")))
  cat("Function value:", x$value, "\n")
  cat("Parameters:\n")
  print(x=unlist(unclass(x$par)), ...)
  cat(sprintf("\n\t%d of %d parameters have been calibrated.\n\n", 
              x$nactive, x$npar))
  cat("Counts:\n")
  print(x$counts)
  cat("Partial fitness values:\n")
  print(x$partial)
  return(invisible())
}

# Methods for optimES.result class ----------------------------------------

#' @export
print.optimES.result = function(x, short=FALSE, ...) {
  
  cat("\nFunction value:", x$value, "\n")
  if(!isTRUE(short)) {
    cat(sprintf("Parameters (%d of %d parameters active).\n",
                length(x$active$par), length(x$par)))
    print(x=x$par, ...)    
    if(!isTRUE(x$active$flag)) cat("* Parameters not calibrated.\n")
  }
  
}

#' @export
coef.optimES.result = function(object, ...) {
  return(object$par)
}

#' @export
plot.optimES.result = function(x, ...) {
  return(invisible(NULL))
}

#' @export
summary.optimES.result = function(object, ...) {
  class(object) = "summary.optimES.result"
  return(object)
}

#' @export 
print.summary.optimES.result = function(x, ...) {
  cat("\nFunction value:", x$value, "\n")
  cat("Parameters:\n")
  print(x=x$par, ...)
  if(!isTRUE(x$active$flag)) cat("* Only active parameters are shown.")

  cat("Partial fitness values:\n")
  print(x$partial)
  
  cat("Counts:\n")
  print(x$counts)
  
}

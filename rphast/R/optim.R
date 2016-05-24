##' Optimize an R function using phast's numerical optimization procedure
##'
##' This function works very much like the optim function in the stats
##' package.  In many phast applications, however, I have noticed that this
##' function converges just as well while taking many fewer function
##' evaluations.  It uses the same optimization routine as phyloFit.  In
##' general it is most efficient to use phyloFit, because some efficiency
##' is lost in passing objects back and forth from R to C (as is necessary
##' when using C code to optimize an R function, whereas phyloFit uses
##' C code to optimize a C function).
##' @title Optimize using phast's optimization code
##' @param func A function to be maximized.  The first argument of the function
##' should be a numeric vector of the parameters to be optimized
##' @param params A vector of initial values to send as the first argument
##' of func.
##' @param lower A vector of the same length as the vector of parameters to
##' be optimized, giving the lower bounds for each parameter.  If NULL,
##' set the lower bounds to -Inf for all parameters.
##' @param upper A vector of the same length as the vector of parameters to
##' be optimized, giving the upper bounds for each parameter.  If NULL,
##' set the upper bounds to Inf for all parameters.
##' @param precision The "precision" to use for the optimization, which
##' affects convergence criteria.  Choices are "LOW", "MED", "HIGH", or
##' "VERY_HIGH".
##' @param logfile If non-NULL, give the name of a file to write an optimization
##' log to
##' @param ... Additional arguments to be passed to func at each function
##' call.  These arguments will not be optimized.
##' @return A list with three elements: value: the optimized value of the
##' function, par: a vector giving the parameters at the optimized value, and
##' neval: the number of function evaluations used in the optimization.
##' @author Melissa J. Hubisz and Adam Siepel
##' @export
optim.rphast <- function(func, params, lower=NULL, upper=NULL, precision="HIGH", logfile=NULL, ...) {
  argl <- list(...)
  if (length(argl) >= 1L) {
    for (i in 1:length(argl)) {
      if (is.msa(argl[[i]])) {
        argl[[i]] <- as.pointer.msa(argl[[i]])
      }
    }
  }
  mylikelihood <- function(x, l=argl) {
    l1 <- list(x)
    if (length(l) >= 1L) {
      for (i in 1:length(l))
        l1[[names(l)[i]]] <- l[[i]]
    }
    names(l1)[1] <- names(formals(func))[1]
    do.call(func, l1)
  }
  rv <- .Call.rphast("rph_opt_bfgs", mylikelihood, params, lower, upper, precision, logfile, new.env())
  rv
}


  

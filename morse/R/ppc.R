#' Posterior predictive check plot
#' 
#' Plot posterior predictive check for \code{reproFitTT} and \code{survFitTT}
#' objects.
#' 
#' @param x an object used to select a method.
#' @param \dots Further arguments to be passed to generic methods.

#' @export
ppc <- function(x, ...){
  UseMethod("ppc")
}

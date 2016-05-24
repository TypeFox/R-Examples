#' Interquartile Range generic
#' 
#' Compute the interquartile range.
#' 
#' @param x an object.
#' @param \dots any additional arguments. These are primarily used in \code{IQR.default}
#'   which calls \code{stats::IQR}.
#' @details If \code{x} is an object of class \code{Bolstad} then the posterior 
#'   IQR of the parameter of interest will be calculated.
#' @author James Curran
#' @export
IQR = function(x, ...){
  UseMethod("IQR")
}

#' @export
IQR.default = function(x, ...){
  stats::IQR(x, ...)
}

#' @export
IQR.Bolstad = function(x, ...){
  return(diff(quantile(x, probs = c(0.25, 0.75), ...)))
}
#' Median generic
#' 
#' @param x an object.
#' @param na.rm Ideally if \code{TRUE} then missing values will be removed, but not currently used
#' @details If \code{x} is an object of class \code{Bolstad} then the posterior 
#'   median of the parameter of interest will be calculated.
#' @author James Curran
#' @method median Bolstad
#' @export
median.Bolstad = function(x, na.rm = FALSE){
  return(quantile(x, probs = 0.5))
}


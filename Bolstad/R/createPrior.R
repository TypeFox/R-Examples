#' Create prior generic
#' 
#' @param x a vector of x values at which the prior is to be specified (the support of the prior).
#' @param \dots optional exta arguments. Not currently used.
#' @return a linear interpolation function where the weights have been scaled so
#'   the function (numerically) integrates to 1.
#' @export
createPrior = function(x, ...){
  UseMethod("createPrior")
}

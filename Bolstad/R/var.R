#' Variance generic
#'
#' @param x an object for which we want to compute the variance
#' @param \dots Any additional arguments to be passed to \code{var}.
#' @export
var = function(x, ...){
  UseMethod("var")
}
#' @export
var.default = function(x, ...){
  stats::var(x, ...)
}
#' @export
var.Bolstad = function(x, ...){
  if(any(grepl("var", names(x))))
    return(x$var)
  
  xVals = x$param.x
  mx = mean(x, ...)
  fx = approxfun(xVals, (xVals - mx)^2 * x$posterior)
  
  return(integrate(fx, min(xVals), max(xVals))$value)
}

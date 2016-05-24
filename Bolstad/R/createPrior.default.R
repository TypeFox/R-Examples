#' Create prior default method
#' 
#' @param x a vector of x values at which the prior is to be specified (the support of the prior).
#' @param wt a vector of weights corresponding to the weight of the prior at the
#'   given x values.
#' @param \dots optional exta arguments. Not currently used.
#' @return a linear interpolation function where the weights have been scaled so
#'   the function (numerically) integrates to 1.
#' @export
#' 
createPrior.default = function(x, wt, ...){
  if(length(x) != length(wt))
    stop("x and wt must be of equal length")
  if(any(wt < 0)){
    stop("All weights must be >= 0")
  }
  fx = approxfun(x, wt, yleft = 0, yright = 0, rule = 2)
  A = integrate(fx, x[1], x[length(x)])$value
  if(A != 1){
    cat("Normalizing prior. Normalizing constant: ")
    cat(A, "\n")
    fx = approxfun(x, wt/A, yleft = 0, yright = 0, rule = 2)
  }
  
  return(fx)
}
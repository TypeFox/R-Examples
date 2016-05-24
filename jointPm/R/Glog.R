#' logistic cdf - tested
#' @export
#' @param x is variable 1
#' @param y is variable 2
#' @param alpha is parameter (0,1)
Glog=function(x,y,alpha) exp(-(x^(-1/alpha) + y^(-1/alpha))^alpha) # The distribution function of logistic model (assumes tranformed to Frechet scale)
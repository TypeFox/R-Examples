#' negative logistic cdf - NOT WORKING
#' @export
#' @param x is variable 1
#' @param y is variable 2
#' @param s is parameter (0,Inf)
Gneglog=function(x,y,alpha) exp(-1/x-1/y+(x^(alpha)+y^(alpha))^(-1/alpha))
 
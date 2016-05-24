#' Numerical integration using Simpson's Rule
#' 
#' Takes a vector of \eqn{x} values and a corresponding set of postive
#' \eqn{f(x)=y} values, or a function, and evaluates the area under the curve:
#' \deqn{ \int{f(x)dx} }.
#' 
#' 
#' @param x a sequence of \eqn{x} values.
#' @param fx the value of the function to be integrated at \eqn{x} or a
#' function
#' @param n.pts the number of points to be used in the integration. If \code{x}
#' contains more than n.pts then n.pts will be set to \code{length(x)}
#' @return A list containing two elements, \code{value} - the value of the
#' intergral, and \code{cdf} - a list containing elements x and y which give a
#' numeric specification of the cdf.
#' @keywords misc
#' @examples
#' 
#' ## integrate the normal density from -3 to 3
#' x = seq(-3, 3, length = 100)
#' fx = dnorm(x)
#' estimate = sintegral(x,fx)$value
#' true.val = diff(pnorm(c(-3,3)))
#' abs.error = abs(estimate-true.val)
#' rel.pct.error =  100*abs(estimate-true.val)/true.val
#' cat(paste("Absolute error :",round(abs.error,7),"\n"))
#' cat(paste("Relative percentage error :",round(rel.pct.error,6),"percent\n"))
#' 
#' ## repeat the example above using dnorm as function
#' x = seq(-3, 3, length = 100)
#' estimate = sintegral(x,dnorm)$value
#' true.val = diff(pnorm(c(-3,3)))
#' abs.error = abs(estimate-true.val)
#' rel.pct.error =  100*abs(estimate-true.val)/true.val
#' cat(paste("Absolute error :",round(abs.error,7),"\n"))
#' cat(paste("Relative percentage error :",round(rel.pct.error,6)," percent\n"))
#' 
#' ## use the cdf
#' 
#' cdf = sintegral(x,dnorm)$cdf
#' plot(cdf, type = 'l', col = "black")
#' lines(x, pnorm(x), col = "red", lty = 2)
#' 
#' ## integrate the function x^2-1 over the range 1-2
#' x = seq(1,2,length = 100)
#' sintegral(x,function(x){x^2-1})$value
#' 
#' ## compare to integrate
#' integrate(function(x){x^2-1},1,2)
#' 
#' 
#' @export sintegral
sintegral = function(x, fx, n.pts = max(256, length(x))){
  ## numerically integrates fx over x using Simpsons rule
  ## x - a sequence of x values
  ## fx - the value of the function to be integrated at x
  ##    - or a function
  ## n.pts - the number of points to be used in the integration
  
  if(class(fx) == "function")
    fx = fx(x)
  
  n.x = length(x)
  
  if(n.x != length(fx))
    stop("Unequal input vector lengths")
  
  ## Shouldn't need this any more
#   if(n.pts < 64)
#     n.pts = 64
  
  ## use linear approximation to get equally spaced x values
  
  
  ap = approx(x, fx, n = 2 * n.pts + 1)
  
  h = diff(ap$x)[1]
  
  integral = h*(ap$y[2 * (1:n.pts) - 1] +
                  4 * ap$y[2 * (1:n.pts)] +
                  ap$y[2 * (1:n.pts) + 1]) / 3
  
  
  results = list(value = sum(integral),
                 cdf = list(x = ap$x[2*(1:n.pts)], y = cumsum(integral)))
  class(results) = "sintegral"
  return(results)
}

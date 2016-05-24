#' Posterior quantiles
#' 
#' @param x an object of class \code{Bolstad}
#' @param probs numeric vector of probabilities with values in \eqn{[0,1]}.
#' @param \dots, any extra arguments needed.
#' @details If \code{x} is of class \code{Bolstad} then this will find the
#' quantiles of the posterior distribution using numerical integration and
#' linear interpolation if necessary.
#' @method quantile Bolstad
#' @export 
quantile.Bolstad = function(x, probs = seq(0, 1, 0.25), ...){
  if(any(grepl("quantileFun", names(x))))
    return(x$quantileFun(probs, ...))
  
  res = sintegral(x$param.x, x$posterior)$cdf
  qFn = approxfun(res$y, res$x)
  
  return(qFn(probs))
}

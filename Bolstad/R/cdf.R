#' Cumulative distribution function generic
#' 
#' This function returns the cumulative distribution function (cdf) of the 
#' posterior distribution of the parameter interest over the range of values for
#' which the posterior is specified.
#' 
#' @param x An object for which we want to compute the cdf
#' @param \dots Any other parameters. Not currently used.
#' @return either the exact cdf of the posterior if a conjugate prior has been
#'   used, or a a \code{stats::splinefun} which will compute the lower
#'   tail probability of the parameter for any valid input.
#' @author James Curran
#' @export 
cdf = function(x, ...){
  UseMethod("cdf")
}

#' @describeIn cdf Cumulative distribution function for posterior density
#' @export 
#' 
cdf.Bolstad = function(x, ...){
  if(class(x) != "Bolstad")
    stop("x must be an object of class Bolstad")
  
  if(any(grepl("cdf", names(x))))
    return(x$cdf)
  
  res = sintegral(x$param.x, x$posterior)$cdf
  return(splinefun(res$x, res$y))
}

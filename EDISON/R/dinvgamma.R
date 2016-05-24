#' Calculate inverse gamma distribution.
#' 
#' This function calculates the density of the inverse gamma distribution.
#' 
#' 
#' @param x Input.
#' @param shape Shape parameter.
#' @param scale Scale parameter (1/rate).
#' @param log Whether to return the log density.
#' @return Returns the density (or log density).
#' @author Frank Dondelinger
#' @examples
#' 
#' # Draw samples from inverse gamma distribution with shape parameter 1 
#' # and scale parameter 1
#' samples = rinvgamma(100, shape=1, scale=1)
#' 
#' # Calculate density of samples
#' densities = dinvgamma(samples, shape=1, scale=1)
#' 
#' @export dinvgamma
dinvgamma <-
function(x, shape, scale=1, log=FALSE) {

    log.dens = shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) 
            - (scale/x)
  
    if(!log) {
      return(exp(log.dens))
    } else {
      return(log.dens)
    }
}


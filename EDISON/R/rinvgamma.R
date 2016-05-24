#' Samples from the inverse gamma distribution.
#' 
#' This function samples from the inverse gamma distribution.
#' 
#' 
#' @param n Number of values to sample.
#' @param shape Shape parameter.
#' @param scale Scale parameter (1/rate).
#' @return Random sample from the inverse gamma distribution.
#' @author Frank Dondelinger
#' @seealso \code{\link{dinvgamma}}
#' @examples
#' 
#' # Draw samples from inverse gamma distribution with shape parameter 1 
#' # and scale parameter 1
#' samples = rinvgamma(100, shape=1, scale=1)
#' 
#' # Calculate density of samples
#' densities = dinvgamma(samples, shape=1, scale=1)
#' 
#' @export rinvgamma
rinvgamma <-
function(n, shape, scale) {
    return(1 / rgamma(n, shape=shape, scale=1/scale))
}


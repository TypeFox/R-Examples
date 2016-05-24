##' Compute mean of truncated normal
##'
##' Given a particular cutoff (\code{x}), what is the mean of the new truncated normal?
##'	
##' The mean of a pdf (f(x)) is computed by integrating the function x*f(x) from -Inf to Inf. To compute
##' a truncated normal, the function remains the same, but the limits change (from -Inf to \code{cutoff}). 
##' @param cutoff The cutoff value for the truncated normal
##' @param mean The mean of the non-truncated distribution (defaults to zero)
##' @param sd The sd of the non-truncated distribution (defaults to 1)
##' @return the mean of the tuncated normal distribution
##' @author Dustin Fife
##' @export
##' @examples
##' ### compute mean of a distribution cutoff at mean
##' normalMean(0, 0, 1)
##' ### compare to simulated data
##' x = rnorm(1000000, 0, 1)
##' mean(x[x>0])
normalMean = function(cutoff, mean=0, sd=1){
	
	### specify function to integrate
	fun = function(x, mean, sd){x*dnorm(x, mean, sd)}
	
	### integrate it
	integral = integrate(f=fun, lower=cutoff, upper=Inf, mean=mean, sd=sd)$value/pnorm(cutoff, lower.tail=F)
	return(integral)	
}
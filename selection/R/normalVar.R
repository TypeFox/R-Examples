##' Compute variance of truncated normal
##'
##' Given a particular cutoff (\code{x}), what is the variance of the new truncated normal?
##'	
##' The variance of a pdf (f(x)) is computed by integrating the function x^2*f(x) from -Inf to Inf. To compute
##' a truncated normal, the function remains the same, but the limits change (from -Inf to \code{cutoff}). 
##' @param cutoff The cutoff value for the truncated normal
##' @param truncMean The mean of the truncated normal. If left NULL, it will be estimated.
##' @param mean The mean of the non-truncated distribution (defaults to zero)
##' @param sd The sd of the non-truncated distribution (defaults to 1)
##' @return the variance of the tuncated normal distribution
##' @author Dustin Fife
##' @export
##' @examples
##' ### compute variance of a distribution cutoff at zero
##' normalVar(1)
##' ### compare to simulated data
##' x = rnorm(10000000, 0, 1)
##' var(x[x>1])
normalVar = function(cutoff, truncMean=NULL, mean=0, sd=1){

	### if they haven't estimated the truncated mean, estimate it for them
	if (is.null(truncMean)){
		truncMean = normalMean(cutoff, mean, sd)
	}
	
	### create function to later integrate
	fun = function(cutoff, truncMean=truncMean, mean, sd){
		((cutoff-truncMean)^2*dnorm(cutoff, mean, sd))
	}
	
	#### integrate then divide to normalize it
	(integrate(fun, cutoff, Inf, truncMean=truncMean, mean=mean, sd=sd)$value)/pnorm(cutoff, lower.tail=F)
}
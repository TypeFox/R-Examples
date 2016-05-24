#' @title 
#' Plot Method for EpiBayes Output
#' 
#' @description
#' This method plots the output of the function \code{\link{EpiBayes_ns}} and 
#'     \code{\link{EpiBayes_s}}. It does so by plotting the posterior distribution(s) of the 
#'     cluster-level prevalence(s).
#'
#' @param x An object of class \code{eb} (e.g., the output of functions 
#'     \code{\link{EpiBayes_ns}} or \code{link{EpiBayes_s}}).
#' @param burnin Number of MCMC iterations to discard from the beginning of the chain. 
#'     Integer scalar.
#' @param ... Additional arguments to be passed to \code{plot}.
#'
#' @return
#' One plotting window for each subzone including posterior distributions for the 
#'     cluster-level prevalence across all time periods.
#'
#' @export
#' @name plot.eb
#' @import scales
plot.eb = function(x, burnin = NULL, ...){
	
	# Extract useful objects from the EpiBayes output
	taumat = x$taumat
	H = x$ForOthers$H
	reps = x$ForOthers$reps
	MCMCreps = x$ForOthers$MCMCreps
	if (is.null(burnin)){ 
		burnin = x$ForOthers$burnin
		}

	# browser()
	
	# Plot the posterior distributions
		# First, get the horizontal axis plotting values
		xdens = density(taumat[1, 1, -c(1:burnin)], from = 0, to = 1)$x
		
	for (i in 1:H){
		
		denstemp = apply(
						matrix(taumat[, i, -c(1:burnin)], reps, MCMCreps - burnin), 
						1, 
						function(z){
									density(z, from = 0, to = 1)$y
									})
									
		plot(
			c(0, 1), 
			c(0, max(denstemp)), 
			type = "n",
			main = "Posterior Distributions for Cluster-level Prevalence \n for each Replicated Data Set",
			xlab = "Cluster-level Prevalence",
			ylab = "Density",
			...
			)
		
		apply(
			denstemp, 
			2, 
			function(z){
						lines(xdens, z, type = "l", col = scales::alpha("black", 0.7))
						})
	}
}
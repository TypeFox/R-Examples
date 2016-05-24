#' @title 
#' Summary Method for EpiBayes Historical Object
#' 
#' @description
#' This function gives summary measurements for posterior distributions of cluster-level 
#'     prevalences across all time periods considered. It does so by examining the object 
#'     output by the \code{\link{EpiBayesHistorical}} function of class
#'     \code{ebhistorical}.
#'
#' @param object An object of class \code{ebhistorical} (e.g., the output of 
#'     function \code{\link{EpiBayesHistorical}}).
#' @param sumstat The summary statistic that the user wishes to be output. The three choices are \code{mean}, \code{quantile}, and \code{variance}. Character scalar.
#' @param prob The probability associated with the quantile one wishes to calculate. Only
#'     used if \code{sumstat} is chosen to be \code{quantile}. Real scalar.
#' @param time.labels Labels for the time period axis (e.g., a vector of years). Character
#'     vector.
#' @param busterapprox Boolean value indicating whether or not the summary statistics 
#'     should be computed using the raw posterior distribution values computed via MCMC 
#'     (\code{TRUE}) or using the best beta distribution approximation via moment matching
#'     (\code{FALSE}). Boolean.
#' @param burnin Number of MCMC iterations to discard from the beginning of the chain. 
#'     Integer scalar.
#' @param ... Additional arguments to be passed on to summary.
#'
#' @return
#' The summary statistics are returned in a matrix with rows indexed by time period and 
#'     columns indexed by subzone. The body of the matrix is filled with the summary 
#'     statistics requested for by the argument \code{sumstat}.
#'
#' @seealso 
#' This is a method for objects of class \code{ebhistorical} returned by the function
#'     \code{\link{EpiBayesHistorical}} and creates its own class of object much like the 
#'     summary method for \code{lm} objects. There is an additional plot method that will 
#'     plot the summary output from this method, \code{summary.ebhistorical}.
#'
#' @export
#' @name summary.ebhistorical
summary.ebhistorical = function(object, sumstat = "mean", prob = 0.95, time.labels = NULL, busterapprox = FALSE, burnin = NULL, ...){
		
		# Extract relevant information
		periodi.tau = object$RawPost
		periodi.betabust = object$BetaBusterEst
		unique.subzones = object$ForOthers$unique.subzones
		if (is.null(burnin)){ 
			burnin = object$ForOthers$burnin
		}
		
		# Get dimensions for output
		n.periods = length(periodi.betabust)
		n.subzones = nrow(periodi.betabust[[1]])
		
		# Compute summary statistic for each subzone
		if(busterapprox){
			sumstat.mat = matrix(unlist(lapply(periodi.betabust,
						   function(x){
						       apply(x, 
						       		 1, 
						       		 function(y){
						       		 	switch(sumstat,
						       		 			"mean" = y[1]/(y[1] + y[2]),
						       		 			"quantile" = qbeta(prob, y[1], y[2]),
						       		 			"variance" = (y[1] * y[2])/((y[1] + y[2])^2 * (y[1] + y[2] + 1))
						       		 	)
						       		 }
						       	)
						    }
						)
					),
					n.subzones, n.periods)
			} else{		
				sumstat.mat = matrix(unlist(lapply(periodi.tau,
							   function(x){
							       apply(x, 
							       		 1, 
							       		 function(y){
							       		 	switch(sumstat,
							       		 			"mean" = mean(y[-c(1:burnin)]),
							       		 			"quantile" = quantile(y[-c(1:burnin)], prob = prob),
							       		 			"variance" = var(y[-c(1:burnin)])
							       		 	)
							       		 }
							       	)
							    }
							)
						),
						n.subzones, n.periods)		
					}	
		
		sumstat.mat = list(
						"sumstat.mat" = sumstat.mat, 
						"sumstat.call" = match.call(expand.dots = TRUE), 
						"ForOthers" = list(
										unique.subzones = unique.subzones,
										time.labels = time.labels
										)
						)
		# Create a summary class for this object so it can be plotted
		class(sumstat.mat) = "ebhistoricalsummary"
		
		# Print output to console
		print(sumstat.mat)
		
		return(invisible(sumstat.mat))
							
	}
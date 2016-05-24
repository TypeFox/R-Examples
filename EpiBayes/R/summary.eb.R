#' @title 
#' Summary Method for EpiBayes Object
#' 
#' @description
#' This function gives summary measurements for posterior distributions of cluster-level 
#'     prevalences across all time periods considered. It does so by examining the object 
#'     output by the \code{\link{EpiBayes_ns}} or \code{\link{EpiBayes_s}} function of 
#'     class \code{eb}.
#'
#' @param object An object of class \code{ebhistorical} (e.g., the output of 
#'     function \code{\link{EpiBayesHistorical}}).
#' @param prob The probability associated with the highest posterior density (HPD) 
#'     intervals one wishes to calculate for each of the reported parameters. 
#' @param burnin Number of MCMC iterations to discard from the beginning of the chain. 
#'     Integer scalar.
#' @param n.output Number of replicated data sets' summary measures to print. Integer 
#'     scalar.
#' @param ... Additional arguments to be passed on to summary.
#'
#' @return
#' The summary statistics are returned in a list with the first entry containing the 
#'     simulation output (\code{p2.tilde}, \code{p4.tilde}, and \code{p6.tilde}), the next 
#'     containing summary measures for the first ten replicated data sets' \code{gam}, and 
#'     the rest containing summary measures for the first ten replicated data sets' 
#'     \code{tau} values (one for each subzone, if applicable).
#'     The summary measurements taken on the posterior distributions include the posterior
#'     mean, standard deviation, standard error of the mean, time-series adjusted standard 
#'     error of the mean, and the lower and upper HPD interval limits, in that order.
#'     For reference purposes, below are the descriptions for the summarized variables.
#' 
#' \tabular{lll}{
#'     Output \tab Description \cr
#'     \code{p2.tilde} \tab Proportion of simulated data sets that result in the probability of \code{poi} prevalence \emph{below} \code{poi.thresh} with probability \code{p1} \cr
#'     \code{p4.tilde} \tab Proportion of simulated data sets that result in the probability of \code{poi} prevalence \emph{above} \code{poi.thresh} with probability \code{p1} \cr
#'     \code{p6.tilde} \tab Proportion of simulated data sets that result in the probability of \code{poi} prevalence \emph{between} \code{poi.lb} and \code{poi.ub} with probability \code{p1} \cr
#'     \code{taumat} \tab Posterior distributions of the cluster-level prevalence for all simulated data sets (i.e., \code{reps}) \cr
#'     \code{gammat} \tab Posterior distribution of the subzone-level prevalence (3-level) OR Posterior distribution of the probability of the disease being in the region (2-level) \cr
#' }
#'
#' @seealso 
#' This is a method for objects of class \code{eb} returned by the function
#'     \code{\link{EpiBayes_ns}} or \code{\link{EpiBayes_s}} and creates its own class of 
#'     object much like the summary method for \code{lm} objects. 
#'
#' @export
#' @name summary.eb
#' @import coda
summary.eb = function(object, prob = 0.95, burnin = NULL, n.output = NULL, ...){
	
		# Extract some values for processing
		reps = object$ForOthers$reps
		MCMCreps = object$ForOthers$MCMCreps
		H = object$ForOthers$H
		p2.tilde = object$p2.tilde
		p4.tilde = object$p4.tilde
		p6.tilde = object$p6.tilde
		if (is.null(burnin)){ 
			burnin = object$ForOthers$burnin
		}
		
		# Set number of output rows to 10 at most (and error checking)
		if (is.null(n.output)){
			if (reps < 10){
				n.output = reps
			} else{
				n.output = 10
			}
		} else if (n.output > reps){
			n.output = 10
		}
		
		# Set up an empty list to store varible-specific output
		epibayessum.out = list()
		
		# Store the simulation output
		epibayessum.out[[1]] = matrix(c(p2.tilde, p4.tilde, p6.tilde), 1, 3)
		colnames(epibayessum.out[[1]]) = c("p2.tilde", "p4.tilde", "p6.tilde")

		# Get summary for gam
		gampost = coda::as.mcmc(t(matrix(object$gammat[, -c(1:burnin)], reps, MCMCreps - burnin)))
		gamsum = matrix(summary(gampost, ...)[[1]], nrow = reps)[1:n.output, ]
		gamHPD = matrix(coda::HPDinterval(gampost, prob = prob)[1:n.output, ], nrow = n.output)
		epibayessum.out[[2]] = matrix(c(gamsum, gamHPD), nrow = n.output)
		colnames(epibayessum.out[[2]]) = c("Mean", "SD", "Naive SE", "Time-series SE", "Lower HPD Limit", "Upper HPD Limit")

		# Get summary for tau
		
		for (i in 1:H){
			taupost.i = coda::as.mcmc(t(matrix(object$taumat[, i, -c(1:burnin)], reps, MCMCreps - burnin)))
			tausum = matrix(summary(taupost.i, ...)[[1]], nrow = reps)[1:n.output, ]
			tauHPD = matrix(coda::HPDinterval(taupost.i, prob = prob)[1: n.output, ], nrow = n.output)
			epibayessum.out[[i + 2]] = matrix(c(tausum, tauHPD), nrow = n.output)
			colnames(epibayessum.out[[i + 2]]) = c("Mean", "SD", "Naive SE", "Time-series SE", "Lower HPD Limit", "Upper HPD Limit")
			}
			
		# Set class of summary object
		class(epibayessum.out) = "ebsummary"
		
		# Print the output to the console
		print(epibayessum.out)
		
		return(invisible(epibayessum.out))
					
	}
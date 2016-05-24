#' @title 
#' Print Method for EpiBayes Object
#' 
#' @description
#' This function prints basic summary measurements for posterior distributions of cluster- 
#'     level prevalences across all time periods considered. It does so by examining the 
#'     object output by the \code{\link{EpiBayes_ns}} or \code{\link{EpiBayes_s}} function 
#'     of class \code{eb}.
#'
#' @param x An object of class \code{eb} (e.g., the output of 
#'     function \code{\link{EpiBayes_ns}}) or \code{\link{EpiBayes_s}}).
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
#'     By default, the function returns summary values for up to ten replicated data sets of
#'     \code{gam} and \code{tau} and reports 95% HPD intervals. 
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
#'     \code{\link{EpiBayes_ns}} or \code{\link{EpiBayes_s}}. The \code{summary} method for
#'     the \code{eb} object class allows for customization of this output.
#'
#' @export
#' @name print.eb
print.eb = function(x, ...){
	
		# Extract some values for processing
		reps = x$ForOthers$reps
		MCMCreps = x$ForOthers$MCMCreps
		H = x$ForOthers$H
		p2.tilde = x$p2.tilde
		p4.tilde = x$p4.tilde
		p6.tilde = x$p6.tilde
		burnin = x$ForOthers$burnin
		
		# Set number of output rows to 10 at most
		if (reps < 10){
			n.output = reps
		} else{
			n.output = 10
		}
		
		# Set up an empty list to store varible-specific output
		epibayessum.out = list()
		
		# Store the simulation output
		epibayessum.out[[1]] = matrix(c(p2.tilde, p4.tilde, p6.tilde), 1, 3)
		colnames(epibayessum.out[[1]]) = c("p2.tilde", "p4.tilde", "p6.tilde")

		# Get summary for gam
		gampost = coda::as.mcmc(t(matrix(x$gammat[, -c(1:burnin)], reps, MCMCreps - burnin)))
		gamsum = matrix(summary(gampost, ...)[[1]], nrow = reps)[1:n.output, ]
		gamHPD = matrix(coda::HPDinterval(gampost, prob = 0.95)[1:n.output, ], nrow = n.output)
		epibayessum.out[[2]] = matrix(c(gamsum, gamHPD), nrow = n.output)
		colnames(epibayessum.out[[2]]) = c("Mean", "SD", "Naive SE", "Time-series SE", "Lower HPD Limit", "Upper HPD Limit")

		# Get summary for tau
		
		for (i in 1:H){
			taupost.i = coda::as.mcmc(t(matrix(x$taumat[, i, -c(1:burnin)], reps, MCMCreps - burnin)))
			tausum = matrix(summary(taupost.i, ...)[[1]], nrow = reps)[1:n.output, ]
			tauHPD = matrix(coda::HPDinterval(taupost.i, prob = 0.95)[1: n.output, ], nrow = n.output)
			epibayessum.out[[i + 2]] = matrix(c(tausum, tauHPD), nrow = n.output)
			colnames(epibayessum.out[[i + 2]]) = c("Mean", "SD", "Naive SE", "Time-series SE", "Lower HPD Limit", "Upper HPD Limit")
			}
					
		print(epibayessum.out)
		
		return(invisible(epibayessum.out))
					
	}
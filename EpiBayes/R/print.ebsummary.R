#' @title 
#' Print Method for EpiBayes Historical Object
#' 
#' @description
#' This function prints the output of objects of class \code{ebhistoricalsummary} from the 
#'     method \code{\link{summary.ebhistorical}} in a nice form.
#'
#' @param x An object of class \code{ebsummary} (e.g., the output of 
#'     function \code{\link{summary.eb}}).
#' @param ... Additional arguments to be passed on to print.
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
#'     Here, \code{poi} refers to the p(arameter) o(f) i(nterest) and is chosen by the user
#'     to be either \code{gam} or \code{tau}. By default, it is chosen to be the cluster-
#'     level prevalence, \code{tau}.
#'
#' @seealso 
#' This is a method for objects of class \code{ebsummary} returned by the method
#'     \code{\link{summary.eb}}, which creates its own class of object much like 
#'     the summary method for \code{lm} objects. 
#'
#' @export
#' @name print.ebsummary
print.ebsummary = function(x, ...){

	# Get number of items in output list
	n.elements = length(x)
	
	# Find if the model is for a 2- or 3-level model
	twolevel.bool = ifelse (n.elements == 3, TRUE, FALSE)
	
	# Print the simulation output
	cat("Simulation output for parameter of interest (poi) \n")
	cat("	*p2.tilde: Percentage of the time the disease is not detected above the disease threshold \n")
	cat("	*p4.tilde: Percentage of the time the disease is detected above the disease threshold \n")
	cat("	*p6.tilde: Percentage of the time the disease is detected between the user-supplied lower and upper bounds of interest \n \n")
	
	rownames(x[[1]]) = ""
	print(x[[1]])
	
	cat("\n ----------------------------------------------------------------------------- \n")
	
	# Print the gam summary output
	if (twolevel.bool){
		cat("gam: Probability of disease being in the region \n")
		} else{
			cat("gam: Subzone-level prevalence \n")
			}
			
	print(x[[2]])
	
	cat("\n ----------------------------------------------------------------------------- \n")
	
	for (i in 1:(n.elements - 2)){
		cat(paste0("tau ", i, ": Cluster-level prevalence in subzone ", i, "\n"))
		print(x[[i + 2]])
		
	cat("\n ----------------------------------------------------------------------------- \n")
		
	}

}
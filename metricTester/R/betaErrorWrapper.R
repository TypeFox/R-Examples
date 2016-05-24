#' Read in and calculate type I and II error rates of a set of beta metric tests
#'
#' Reads in the results of a set of betaMultiLinker results and summarizes the type I and
#' II error rates of the various metric-null combinations under different assembly
#' processes.
#'
#' @param working.directory Optional character string specifying the working directory.
#' If missing, the current working directory will be used. 
#' 
#' @details The internal Reduce call in this function is not formatted to my liking (ETM).
#' Ideally, it would return a list of the length of the number of spatial simulations,
#' then a data frame for each spatial simulation. Instead, consolidates all to a single
#' data frame. It would be easy to make it work with for loops, but it should also be
#' possible to revise the Reduce call and make it return results properly formatted. 
#'
#' @return A data frame of spatial simulations, null models, and metric combinations,
#' summarizing the results of all the betaMultiLinker runs in the working directory.
#' Numbers refer to total number of runs that resulted in an error. Errors for the random
#' spatial simulation refer to type I errors, errors for habitat filtering simulation
#' refer to type II errors for detecting filtering, and errors for competitive
#' exclusion refer to type II errors for detecting competition.
#'
#' @export
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.

betaErrorWrapper <- function(working.directory)
{
	temp <- readIn(working.directory)
	checkedList <- lapply(temp, betaErrorChecker)
	temp2 <- Reduce(function(y, z) mapply(rbind, y, z, SIMPLIFY=F),
		lapply(checkedList, unlist, recursive=F), accumulate=F)
	results <- lapply(seq_along(temp2), function(x) apply(temp2[[x]], 2, sum))
	names(results) <- names(temp2)
	results <- t(as.data.frame(results))
	results
}

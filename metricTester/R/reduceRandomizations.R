#' Reduce randomized results to a manageable list of dataframes
#'
#' The metricsNnulls function creates lists of lists of dataframes. This function will
#' combine the dataframes from each null model into a single data frame. The output is a
#' more manageable list of dataframes. 
#'
#' @param randomizations.list The results of a call to metricsNnulls()
#'
#' @details Given a list of lists of dataframes, such as those that come from a call to
#' metricsNnulls, where the first level of the list relates to a given randomization, and
#' each second level is a data frame containing the calculated metrics after randomization
#' according to a given null model, reduces the results to a simpler list of data frames,
#' where each data frame contains all the results from a given null model from the input
#' randomizations.list.
#'
#' @export
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' #simulate tree with birth-death process
#' tree <- geiger::sim.bdtree(b=0.1, d=0, stop="taxa", n=50)
#'
#' sim.abundances <- round(rlnorm(5000, meanlog=2, sdlog=1)) + 1
#'
#' cdm <- simulateComm(tree, richness.vector=10:25, abundances=sim.abundances)
#'
#' #below not run for timing issues on CRAN
#' #rawResults <- metricsNnulls(tree, cdm)
#'
#' #results <- reduceRandomizations(rawResults)

reduceRandomizations <- function(randomizations.list)
{
	#this command successively combines each element from the long list together via an
	#inner anonymous function that mapply(rbinds) things
	finalResults <- Reduce(function(y, z) mapply(rbind, y, z, SIMPLIFY=FALSE), 
		randomizations.list, accumulate=FALSE)
	finalResults
}

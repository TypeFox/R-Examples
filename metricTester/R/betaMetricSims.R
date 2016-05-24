#' Calculate beta metrics under specified tree and community parameters
#'
#' Takes a specified set of tree size, shape, community richness and abundance parameters,
#' and calculates the beta-level phylogenetic community structure metrics.
#'
#' @param tree.size Number of species desired in the total tree.
#' @param richness.vector Number of species to be placed in each plot. See details.
#' @param delta A value for the delta transformation (Pagel 1999). Values greater than 1
#' push the branching events towards the root, while values less than 1 push the branching
#' events closer to the tips. See details for particularly low delta values.
#' @param abundances Vector of abundances, e.g. a repeated series of 1s for a presence/absence
#' community data matrix, a log-normal distribution, etc. See examples.
#' @param beta.iterations Because the type of beta-level phylogenetic community structure
#' metrics used here return a single value per community data matrix, it is not possible
#' to look for inter-metric correlations with only a single matrix and tree. To deal with
#' this, the same tree can be used with different community data matrices. This argument
#' specifies the number of matrices to be used per tree.
#' 
#' @details The richness.vector (number of species to be placed into each plot) is
#' flexible. For instance, one might want give it 10:19, which would create 10 plots
#' of species richness 10, 11, ... 19. But one could also provide rep(10, 10) to create 10
#' plots of 10 species each. If given a small value, e.g. 0.1, the delta parameter
#' (tree shape) can occasionally result in oddly formatted trees that would cause errors.
#' To deal with this, there is an internal check that will recreate a new tree and
#' re-scale it with the desired delta. This has not been tested at delta < 0.1, and is
#' currently programmed with a while loop. Care should be taken not to get R stuck in an
#' indefinite loop at delta values even lower than 0.1
#'
#' @return A data frame of calculated alpha metrics, and the associated species richness
#' of each plot.
#'
#' @export
#'
#' @importFrom ape is.ultrametric
#' @importFrom geiger sim.bdtree rescale
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' test <- betaMetricSims(tree.size=50, richness.vector=30:40, delta=1,
#'	abundances=round(rlnorm(5000, meanlog=2, sdlog=1)) + 1, beta.iterations=10)

betaMetricSims <- function(tree.size, richness.vector, delta, abundances, beta.iterations)
{
	#simulate tree with birth-death process
	tree <- sim.bdtree(b=0.1, d=0, stop="taxa", n=tree.size)

	newTree <- rescale(tree, "delta", delta)

	ok <- ape::is.ultrametric(newTree)

	#was having trouble with some trees and delta parameters making trees that were not
	#ultrametric and would throw errors when trying to use. this will make a new tree
	#if the first is not ultrametric
	while(ok==FALSE)
	{
		tree <- sim.bdtree(b=0.1, d=0, stop="taxa", n=tree.size)
		newTree <- rescale(tree, "delta", delta)
		ok <- ape::is.ultrametric(newTree)
	}

	#results <- list()
	results <- matrix(nrow=beta.iterations, ncol=length(defineBetaMetrics()))
	
	for(i in 1:beta.iterations)
	{
		cdm <- simulateComm(newTree, richness.vector=richness.vector,
			abundances=abundances)
		prepped <- prepData(newTree, cdm)
		results[i,] <- as.matrix(calcBetaMetrics(prepped)[1,])
	}
	
	colnames(results) <- names(defineBetaMetrics())
	rownames(results) <- paste("beta.iteration", 1:beta.iterations, sep="")
	
	results <- as.data.frame(results)
	
	results
}

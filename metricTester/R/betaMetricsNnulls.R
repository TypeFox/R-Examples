#' Parallelized function that calculates beta metrics on randomized matrices
#'
#' This function sends out jobs to as many cores as are specified. Each randomizes the
#' input CDM according to all defined null models, then calculates each observed beta
#' metric on each randomized matrix.
#'
#' @param tree Phylo object
#' @param picante.cdm A picante-style community data matrix with sites as rows, and
#' species as columns
#' @param optional.dists A symmetric distance matrix can be directly supplied. This option
#' is experimental. Behavior depends on metric being used. If the metric in question
#' relies on the phylogenetic distance matrix from a call to cophenetic(tree), then this 
#' optional distance matrix will be inserted instead.
#' @param regional.abundance A character vector in the form "s1, s1, s1, s2, s2, s3, etc".
#' Optional, will be generated from the input CDM if not provided.
#' @param distances.among A symmetric distance matrix, summarizing the distances among all
#' plots from the cdm. Optional, only used by some null models.
#' @param randomizations The number of times the input CDM should be randomized and the
#' metrics calculated across it.
#' @param cores This function can run in parallel. In order to do so, the user must
#' specify the desired number of cores to utilize. The default is "seq", which runs the
#' calculations sequentially.
#' @param nulls Optional list of named null model functions to use. If invoked, this 
#' option will likely be used to run a subset of the defined null models.
#' @param metrics Optional list of named metric functions to use. If invoked, this option
#' will likely be used to run a subset of the defined metrics.
#'
#' @details This function sends out jobs to as many cores as are specified. Each 
#' randomizes the input CDM according to all defined null models, then calculates each
#' observed metric on each randomized matrix.
#'
#' @return A list of lists of vectors. The first level has as many elements as there
#' are randomizations. The second level has one list for each null model. Each element of
#' this second level is a named vector corresponding to the calculated metric at each
#' plot.
#'
#' @export
#'
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doParallel registerDoParallel
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
#' rawResults <- betaMetricsNnulls(tree, cdm, randomizations=3,
#'	nulls=list("richness"=metricTester:::my_richnessNull,
#'	"frequency"=metricTester:::my_frequency))

betaMetricsNnulls <- function(tree, picante.cdm, optional.dists=NULL,
	regional.abundance=NULL, distances.among=NULL, randomizations=2, cores="seq", nulls,
	metrics)
{
	#if a list of named metric functions is not passed in, assign metrics to be NULL, in
	#which case all metrics will be calculated
	if(missing(metrics))
	{
		metrics <- NULL
	}
		
	#if a list of named nulls functions is not passed in, assign nulls to be NULL, in
	#which case all nulls will be run
	if(missing(nulls))
	{
		nulls <- NULL
	}	

	if(cores == "seq")
	{
		#warn that the analysis is being run sequentially
		warning("Not running analysis in parallel. See 'cores' argument.", call.=FALSE)

		#prep the inputs for parallel randomizations
		nullsPrepped <- prepNulls(tree, picante.cdm, regional.abundance, distances.among)
		#call the parallel for loop. each iteration, save a new list of lists, where each
		#inner element are the metrics for a given null model
		randomResults <- foreach(i = 1:randomizations) %do%
		{
			#run the nulls across the prepped data. this randomizes the CDMs all at once
			randomMatrices <- runNulls(nullsPrepped, nulls)
			#prep the randomized CDMs to calculate the metrics across them
			randomPrepped <- lapply(randomMatrices, function(x) 
				prepData(tree=tree, picante.cdm=x, optional.dists=optional.dists))
			#calculate the metrics
			lapply(randomPrepped, calcBetaMetrics, metrics)
		}
	}
	
	if(cores != "seq")
	{
		registerDoParallel(cores)

		#prep the inputs for parallel randomizations
		nullsPrepped <- prepNulls(tree, picante.cdm, regional.abundance, distances.among)
		#call the parallel for loop. each iteration, save a new list of lists, where each
		#inner element are the metrics for a given null model
		randomResults <- foreach(i = 1:randomizations) %dopar%
		{
			#run the nulls across the prepped data. this randomizes the CDMs all at once
			randomMatrices <- runNulls(nullsPrepped, nulls)
			#prep the randomized CDMs to calculate the metrics across them
			randomPrepped <- lapply(randomMatrices, function(x) 
				prepData(tree=tree, picante.cdm=x, optional.dists=optional.dists))
			#calculate the metrics
			lapply(randomPrepped, calcBetaMetrics, metrics)
		}
		registerDoSEQ()
	}

	randomResults
}

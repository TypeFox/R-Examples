#' Calculate phylogenetic community structure beta metrics
#'
#' Given a prepped metrics.input object, calculate all phylogenetic community structure
#' metrics of interest.
#'
#' @param metrics.input Prepped metrics.input object
#' @param metrics Optional list of named metric functions to use. If invoked, this option
#' will likely be used to run a subset of the defined metrics.
#' 
#' @details This function first confirms that the input is of class metrics.input and,
#' if so, then confirms that the metrics to be calculated are in a named list (via
#' checkMetrics), then lapplies all metric functions to the input metrics.input object.
#'
#' @return A data frame with the calculated metrics and the associated species richness
#' and total abundance of all input "communities".
#'
#' @export
#'
#' @importFrom spacodiR spacodi.calc
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' #simulate tree with birth-death process
#' tree <- geiger::sim.bdtree(b=0.1, d=0, stop="taxa", n=50)
#'
#' sim.abundances <- round(rlnorm(5000, meanlog=2, sdlog=1))
#'
#' cdm <- simulateComm(tree, richness.vector=10:25, abundances=sim.abundances)
#'
#' prepped <- prepData(tree, cdm)
#'
#' results <- calcBetaMetrics(prepped)
#'
#' #an example of how to define ones own metrics for use in the metricTester framework
#' #this "metric" simply calculates the richness of each plot in the CDM
#' exampleMetric <- function(metrics.input)
#' {
#'	output <- mean(apply(metrics.input$picante.cdm, 1, lengthNonZeros))
#'	output
#' }
#'
#' calcBetaMetrics(prepped, metrics=list("example"=exampleMetric))

calcBetaMetrics <- function(metrics.input, metrics)
{
	if(!inherits(metrics.input, "metrics.input"))
	{
		stop("Input needs to be of class 'metrics.input'")
	}
	
	#if a list of named metric functions is not passed in, assign metrics to be NULL, in
	#which case all metrics will be calculated
	if(missing(metrics))
	{
		metrics <- NULL
	}
		
	metrics <- checkBetaMetrics(metrics)
	
	tempResults <- lapply(metrics, function(x) x(metrics.input))
	
	#convert the list to a data frame
	results <- as.data.frame(tempResults)

	results
}

#' Generate expectations for null+metric combinations
#'
#' Will generate mean, SD and CI expectations for the desired metric + null combinations.
#'
#' @param picante.cdm A picante-style community data matrix with sites as rows, and
#' species as columns
#' @param tree Phylo object
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
#' specify the desired number of cores to utilize.
#' @param metrics Optional list of named metric functions to use. If invoked, this option
#' will likely be used to run a subset of the defined metrics.
#' @param nulls Optional list of named null model functions to use. If invoked, this 
#' option will likely be used to run a subset of the defined null models.
#' @param concat.by Whether to concatenate the randomizations by richness, plot or both
#' @param output.raw Default is FALSE. Set to TRUE if raw randomized values are preferred
#' (as opposed to summarized mean, SD, CI, etc).
#'
#' @details Given a list of desired metrics (which should always include richness) and
#' null models, will generate the expected mean, standard deviation and confidence
#' intervals based on the number of specified randomizations. This function is flexible in
#' that new metrics and nulls can be added and tested with it. By setting output.raw to
#' TRUE, the function can also output raw randomized values as opposed to the summarized
#' values.
#'
#' @return A list of data frames, where each data frame corresponds to a given null model,
#' and contains the mean, SD, and CI for each metric for that null model.
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
#' #simulate a log normal abundance distribution
#' sim.abundances <- round(rlnorm(5000, meanlog=2, sdlog=1)) + 1
#'
#' #simulate a community of varying richness
#' cdm <- simulateComm(tree, richness.vector=10:13, abundances=sim.abundances)
#'
#' #below not run for timing issues on CRAN
#' #test <- expectations(picante.cdm=cdm, tree=tree, optional.dists=NULL,
#'	#regional.abundance=NULL, distances.among=NULL, randomizations=3, cores="seq",
#'	#nulls=list("richness"=metricTester:::my_richnessNull), 
#'	#metrics=list("richness"=metricTester:::my_richness, "NAW_MPD"=metricTester:::naw_mpd),
#'	#concat.by="both", output.raw=FALSE)
#'
#' #an example of how to explore behavior of a new metric in the metricTester framework
#' #this "metric" simply calculates the richness of each plot in the CDM
#' exampleMetric <- function(metrics.input)
#' {
#'	output <- apply(metrics.input$picante.cdm, 1, lengthNonZeros)
#'	output
#' }
#'
#' #test2 <- expectations(picante.cdm=cdm, tree=tree, optional.dists=NULL,
#'	#regional.abundance=NULL, distances.among=NULL, randomizations=3, cores=1,
#'	#nulls=list("frequency"=metricTester:::my_frequency), 
#'	#metrics=list("richness"=metricTester:::my_richness, "exampleMetric"=exampleMetric),
#'	#concat.by="both", output.raw=FALSE)

expectations <- function(picante.cdm, tree, optional.dists=NULL, regional.abundance=NULL, 
	distances.among=NULL, randomizations, cores, metrics, nulls, concat.by,
	output.raw=FALSE)
{
	#calculate the raw randomized results across whatever metrics and nulls are called
	rawResults <- metricsNnulls(picante.cdm=picante.cdm, tree=tree, 
		optional.dists=optional.dists, regional.abundance=regional.abundance, 
		distances.among=distances.among, randomizations=randomizations, cores=cores,
		metrics=metrics, nulls=nulls)
	#reduce these randomizations into more managable results
	results <- reduceRandomizations(rawResults)
	#if the user wants the raw data, provide that
	if(output.raw==TRUE)
	{
		return(results)
	}
	#otherwise return the summarized data
	else
	{
		output <- lapply(results, summaries, concat.by=concat.by)
		return(output)
	}
}

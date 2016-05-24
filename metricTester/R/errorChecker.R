#' Wrapper for summarizing randomizations and testing significance of observed metrics
#'
#' Given a data frame of observed metrics and a list of randomizations based on different
#' null models, returns a list of data frames summarizing the significance of observed 
#' metrics both at the single plot and the entire arena level.
#'
#' @param observed Data frame of observed metric scores, such as from observedMetrics()
#' @param reduced.randomizations List of random, reduced results, such as those from
#' reduceRandomizations()
#' @param concat.by Whether to concatenate the randomizations by richness, plot or both
#' @param metrics Optional list of named metric functions to use. If invoked, this option
#' will likely be used to run a subset of the defined metrics.
#' 
#' @details This function wraps a number of smaller functions into a useful type I and II
#' error checker. It takes a reduced list of randomizations such as those reduced from
#' metricsNnulls with reduceRandomizations, summarizes the mean,
#' SD, and CI of each metric plus null model either at the richness or plot level,
#' then compares the observed metric scores to those summarized metrics. It return a list
#' with two elements. The first is a list of data frames, where each corresponds to the 
#' standardized effect scores of the observed metrics for a given null model. The second
#' is a list of data frames, where each corresponds to whether a given plot deviates
#' beyond CI. For the latter, 0 corresponds to within CI, 1 corresponds to less than the
#' CI, and 2 corresponds to greater than the CI.
#'
#' @return A list of lists of data frames.
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
#' #run the metrics and nulls combo function
#' #rawResults <- metricsNnulls(tree, cdm, randomizations=3)
#'
#' #summarize the results
#' #results <- reduceRandomizations(rawResults)
#'
#' #calculate the observed metrics from the input CDM
#' #observed <- observedMetrics(tree, cdm)
#'
#' #not run
#' #test <- errorChecker(observed, results, "richness")

errorChecker <- function(observed, reduced.randomizations, concat.by, metrics)
{
	#if a list of named metric functions is not passed in, assign metrics to be NULL, in
	#which case all length of all metrics will be used
	if(missing(metrics))
	{
		metrics <- NULL
	}

	#lapply the summaries function over the reduced randomizations
	summarized <- lapply(reduced.randomizations, summaries, concat.by)
	
	#this is an important command. depending on what you concatenated by, there should be
	#only a single matching column betweeen the tables (either richness or plot), and
	#so it should merge on that. there are not currently any checks for missing values
	#e.g. with regional null, so need to build something in soon
	if(concat.by == "richness" | concat.by== "plot")
	{
		merged <- lapply(summarized, merge, observed)
	}
	else if(concat.by == "both")
	{
		toFeed <- names(summarized)
		merged <- lapply(1:length(toFeed), function(x) 
			list("richness"=merge(summarized[[x]]$richness, observed),
			"plot"=merge(summarized[[x]]$plot, observed)))
		names(merged) <- toFeed
	}

	#this will return a list of data frames, one for each null model, where the first col
	#is whatever we summarized on, and each successive column is the SES of the observed
	#score based on the randomizations
	if(concat.by == "richness" | concat.by== "plot")
	{
		sesResults <- lapply(1:length(merged), function(x) arenaTest(merged[[x]],
			concat.by, metrics))
	}
	
	#this will return a list of lists of data frames. each of first element of lists
	#corresponds to a null model. then within that there is one data frame for richness
	#and one for plot
	else if(concat.by == "both")
	{
		sesResults <- lapply(1:length(toFeed), function(x) 
			list("richness"=arenaTest(merged[[x]]$richness, concat.by="richness",
			metrics), "plot"=arenaTest(merged[[x]]$plot, concat.by="plot", 
			metrics)))
	}
	#set the names right (works for either both or rich/quad)
	names(sesResults) <- names(merged)
	
	#this will return a list of data frames, one for each null model, where the first col
	#is whatever we summarized on, and each successive column is an indicator of whether
	#the observed score in a plot was bigger or lesser
	if(concat.by == "richness" | concat.by== "plot")
	{
		plotResults <- lapply(1:length(merged), function(x) plotTest(merged[[x]], 
			concat.by, metrics))
	}
	
	#structure follows sesResults with argument concat.by=both above
	else if(concat.by == "both")
	{
		plotResults <- lapply(1:length(toFeed), function(x) 
			list("richness"=plotTest(merged[[x]]$richness, concat.by="richness",
			metrics), "plot"=plotTest(merged[[x]]$plot, concat.by="plot", 
			metrics)))
	}

	#set the names right (works for either both or rich/quad)
	names(plotResults) <- names(merged)
	
	results <- list("ses"=sesResults, "plot"=plotResults)
	results
}

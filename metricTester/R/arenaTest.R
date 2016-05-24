#' Calculate SES of each observed metric + null model combination
#'
#' Given a table of results, where means, SDs, and CIs are bound to the observed scores at
#' the corresponding richness or plot, this function calculates standardized effect
#' scores for each observed metric + null model combination. This is intended to be used
#' to test whether observed values deviate beyond expectations based on the distribution
#' of SES per arena.
#'
#' @param results.table Data frame of observed metrics with expected mean, SD and CI bound
#' in. See example.
#' @param concat.by Whether to concatenate results by richness, plot or both. If
#' richness, observed scores are compared to all randomized scores where the plot had
#' the corresponding richness. If plot, observed scores (e.g. those from plot 1)
#' are compared to all randomized plot 1 scores. If both, both are run and each is
#' saved as a separate data frame in a single list.
#' @param metrics Optional list of named metric functions to use. If invoked, this option
#' will likely be used to run a subset of the defined metrics.
#' 
#' @details Given a table of results, where means, SDs, and CIs are bound to the observed
#' scores at the corresponding richness or plot, this function calculates standardized
#' effect scores for each observed metric + null model combination. 
#'
#' @return A data frame of standardized effect scores.
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
#' cdm <- simulateComm(tree, richness.vector=10:25, abundances=sim.abundances)
#'
#' #below not run for example timing issues on CRAN
#'
#' #run the metrics and nulls combo function
#' #rawResults <- metricsNnulls(tree=tree, picante.cdm=cdm, randomizations=2, cores="seq")
#'
#' #reduce the randomizations to a more manageable format
#' #reduced <- reduceRandomizations(rawResults)
#'
#' #calculate the observed metrics from the input CDM
#' #observed <- observedMetrics(tree, cdm)
#'
#' #summarize the means, SD and CI of the randomizations
#' #summarized <- lapply(reduced, summaries, concat.by="richness")
#'
#' #merge the observations and the summarized randomizations to facilitate significance
#' #testing
#' #merged <- lapply(summarized, merge, observed)
#'
#' #calculate the standardized scores of each observed metric as compared to the richness
#' #null model randomization.
#' #arenaTest(merged$richness, "richness")
#'
#' #do the same as above but across all null models. not run
#' #temp <- lapply(1:length(merged), function(x) arenaTest(merged[[x]], "richness"))

arenaTest <- function(results.table, concat.by, metrics)
{
	#if a list of named metric functions is not passed in, assign metrics to be NULL, in
	#which case all length of all metrics will be used
	if(missing(metrics))
	{
		metrics <- NULL
	}
		
	#take advantage of the checkMetrics function to find the name of all metrics
	#get rid of the "richness" metric name
	metricNames <- names(checkMetrics(x=metrics))[2:length(names(checkMetrics(x=metrics)))]
	ses <- list()
	for(i in 1:length(metricNames))
	{
		average.name <- paste(metricNames[i], "average", sep=".")
		SD.name <- paste(metricNames[i], "sd", sep=".")
		average <- results.table[,average.name]
		observed <- results.table[,metricNames[i]]
		#it is possible to have a SD of 0, particularly if a given richness is
		#infrequently sampled in the randomizations, or if only a few randomizations are
		#run. if that happens, clearly the SES would be NA. set the SDs that are 0 to the
		#mean of the SDs
		SD <- results.table[,SD.name]
		SD[SD == 0] <- mean(SD)
		ses[[i]] <- (observed-average)/SD
	}
	if(concat.by=="richness")
	{
		ses <- as.data.frame(ses)
		names(ses) <- metricNames
		ses <- data.frame(richness=results.table$richness, ses)
	}
	else if(concat.by=="plot")
	{
		ses <- as.data.frame(ses)
		names(ses) <- metricNames
		ses <- data.frame(plot=results.table$plot, ses)
	}
	ses
}

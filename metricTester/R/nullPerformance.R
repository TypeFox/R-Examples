#' Summarize null model performance of a series of summarized simulation results
#'
#' Flexible function that summarizes null performance after reading in and testing
#' simulation results with either sesIndiv or plotOverall.
#'
#' @param summarized.results The results of a call to sesIndiv() or plotOverall()
#' @param simulations Default is "all". Alternatively, can supply a vector of simulation
#' names to summarize the results over.
#' @param metrics Default is "all". Alternatively, can supply a vector of metric
#' names to summarize the results over.
#' @param concat.by Default is "both". Alternatively, can supply either "plot" or
#' "richness".
#' 
#' @details If an overall picture of null performance is desired, this function can
#' provide it. It can also be used to summarize null performance over a specific subset
#' of simulations, metrics, and concatenation options. If provided with the results
#' of a call to plotOverall, the options are more limited. Currently, if provided with
#' such a result, the assumption is
#' that there are three spatial simulations, "random", "filtering", and "competition". It
#' then assumes that any clustered or overdispersed plots for the random simulation,
#' or any overdispersed or clustered for the filtering or competition simulations,
#' respectively, count as typeI errors. It assumes that any plots that are not
#' clustered or overdispersed for the filtering or competition simulations, respectively,
#' count as typeII errors.
#'
#' @return A data frame of summarized results
#'
#' @export
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' #not run
#' #results <- readIn()
#' #summ <- sesIndiv(results)
#' #examp <- nullPerformance(summ)

nullPerformance <- function(summarized.results, simulations="all", metrics="all",
	concat.by="both")
{
	if(simulations=="all")
	{
		simulations <- unique(summarized.results$simulation)
	}
	#complicated if statement here. it says if the specified simulations do not contain
	#an entry of character "all" and there is any difference between the specified sims
	#and the unique simulations in the summarized.results table, throw an error
	else if(all(!(simulations %in% "all")) &
		length(setdiff(simulations, unique(summarized.results$simulation))) > 0)
	{
		stop("Specified simulations do not match those in the results table")
	}
	else
	{
		simulations <- simulations
	}

	if(concat.by=="both")
	{
		concat.by <- unique(summarized.results$concat.by)
	}
	else if(all(!(concat.by %in% "both")) &
		length(setdiff(concat.by, unique(summarized.results$concat.by))) > 0)
	{
		stop("concat.by must be set to plot, richness, or both")
	}
	else
	{
		concat.by <- concat.by
	}

	if(metrics=="all")
	{
		metrics <- unique(summarized.results$metric)
	}
	else if(all(!(metrics %in% "all")) &
		length(setdiff(metrics, unique(summarized.results$metric))) > 0)
	{
		stop("Specified metrics do not match those in the results table")
	}
	else
	{
		metrics <- metrics
	}

	nulls <- unique(summarized.results$null.model)
	typeI <- c()
	typeII <- c()

	#use these for sesIndiv results	
	if(names(summarized.results)[5] == "total.runs")
	{
		for(i in 1:length(nulls))
		{
			temp <- summarized.results[summarized.results$null.model %in% nulls[i]
				& summarized.results$simulation %in% simulations
				& summarized.results$concat.by %in% concat.by
				& summarized.results$metric %in% metrics,]
			temp$typeIrate <- 100 * temp$typeI/temp$total.runs
			temp$typeIIrate <- 100 * temp$typeII/temp$total.runs
			typeI[i] <- mean(temp$typeIrate, na.rm=TRUE)
			typeII[i] <- mean(temp$typeIIrate, na.rm=TRUE)
		}
	}

	#use this for plotOverall results
	else if(names(summarized.results)[5] == "clustered")
	{
		#generate some simulation-specific data frames
		random <- summarized.results[summarized.results$metric %in% metrics
			& summarized.results$simulation == "random"
			& summarized.results$concat.by %in% concat.by
			& summarized.results$null.model %in% nulls,]
		filtering <- summarized.results[summarized.results$metric %in% metrics
			& summarized.results$simulation == "filtering"
			& summarized.results$concat.by %in% concat.by
			& summarized.results$null.model %in% nulls,]
		competition <- summarized.results[summarized.results$metric %in% metrics
			& summarized.results$simulation == "competition"
			& summarized.results$concat.by %in% concat.by
			& summarized.results$null.model %in% nulls,]

		#define typeI & II error rates for each of these
		random$typeIrate <- 100 *
			(random$clustered + random$overdispersed)/random$total.plots
		random$typeIIrate <- NA

		filtering$typeIrate <- 100 * filtering$overdispersed/filtering$total.plots
		#notice this important revision. the typeII rate is not simply the number of
		#plots that were not clustered. if the plot was overdispersed than it
		#already counted as a typeI error. revising here and below to incorporate that 
		filtering$typeIIrate <- 100 * 
			(filtering$total.plots - filtering$clustered -
			filtering$overdispersed)/filtering$total.plots

		competition$typeIrate <- 100 *
			competition$clustered/competition$total.plots
		competition$typeIIrate <- 100 * 
			(competition$total.plots -
			competition$overdispersed - competition$clustered)/competition$total.plots

		#redefine summarized.results
		summarized.results <- rbind(random, filtering, competition)

		for(i in 1:length(nulls))
		{
			temp <- summarized.results[summarized.results$null.model %in% nulls[i]
				& summarized.results$simulation %in% simulations
				& summarized.results$concat.by %in% concat.by
				& summarized.results$metric %in% metrics,]
			typeI[i] <- mean(temp$typeIrate, na.rm=TRUE)
			typeII[i] <- mean(temp$typeIIrate, na.rm=TRUE)
		}
	}

	results <- data.frame(nulls, typeI, typeII)
	results
}

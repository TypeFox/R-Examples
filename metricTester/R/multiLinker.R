#' Run multiple simulations and calculations to test metric + null performance
#'
#' This function runs multiple iterations of the linker function, saving results to file.
#'
#' @param no.taxa The desired number of species in the input phylogeny
#' @param arena.length A numeric, specifying the length of a single side of the arena
#' @param mean.log.individuals Mean log of abundance vector from which species abundances
#' will be drawn
#' @param length.parameter Length of vector from which species' locations are drawn. Large
#' values of this parameter dramatically decrease the speed of the function but result in
#' nicer looking communities
#' @param sd.parameter Standard deviation of vector from which species' locations are 
#' drawn
#' @param max.distance The geographic distance within which neighboring
#' indivduals should be considered to influence the individual in question
#' @param proportion.killed The percent of individuals in the total arena that should be
#' considered (as a proportion, e.g. 0.5 = half)
#' @param competition.iterations Number of generations over which to run competition 
#' simulations
#' @param no.plots Number of plots to place
#' @param plot.length Length of one side of desired plot
#' @param concat.by Whether to concatenate the randomizations by richness, plot or both
#' @param randomizations The number of randomized CDMs, per null, to generate. These are
#' used to compare the significance of the observed metric scores.
#' @param cores The number of cores to be used for parallel processing.
#' @param iterations The number of complete tests to be run. For instance, 1 iteration
#' would be considered a complete cycle of running all spatial simulations, randomly
#' placing plots in the arenas, sampling the contents, creating a community data
#' matrix, calculating observed metric scores, then comparing these to the specified
#' number of randomizations of the original CDMs. 
#' @param prefix Optional character vector to affix to the output RData file names, e.g.
#' "test1". 
#' @param simulations Optional list of named spatial simulation functions to use. These
#' must be defined in the defineSimulations function. If invoked, this option will likely
#' be used to run a subset of the defined spatial simulations.
#' @param nulls Optional list of named null model functions to use. If invoked, this 
#' option will likely be used to run a subset of the defined null models.
#' @param metrics Optional list of named metric functions to use. If invoked, this option
#' will likely be used to run a subset of the defined metrics.
#' 
#' @details This function wraps a number of other wrapper functions into
#' one big metric + null performance tester function. Unlike the basic linker function,
#' multiple tests can be run, with results saved as RDS files.
#'
#' @return A list of lists of data frames. The first level of the output has one element 
#' for each simulation. The second level has one element for each null model. Each of
#' these elements is a list of two data frames, one that summarizes the plot-level
#' significance and another and arena-level significance.
#'
#' @export
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' #not run
#' #system.time(multiLinker(no.taxa=50, arena.length=300, mean.log.individuals=3.2, 
#' 	#length.parameter=5000, sd.parameter=50, max.distance=20, proportion.killed=0.3, 
#'	#competition.iterations=2, no.plots=20, plot.length=30, concat.by="richness", 
#'	#randomizations=3, cores="seq", iterations=2, prefix="test",
#'	#nulls=list("richness"=metricTester:::my_richnessNull,
#'	#"frequency"=metricTester:::my_frequency)))

multiLinker <- function(no.taxa, arena.length, mean.log.individuals, length.parameter, 
	sd.parameter, max.distance, proportion.killed, competition.iterations, no.plots, 
	plot.length, concat.by, randomizations, cores, iterations, prefix,
	simulations, nulls, metrics)
{
	#create a simple file name, specify whether concatenating by plot or richness
	for(i in 1:iterations)
	{
		if(is.null(prefix))
		{
			filename <- paste("iteration", i, "by_", concat.by, ".RDS", sep="")
		}
		else
		{
			filename <- paste(prefix, "_", "iteration", i, "by_", concat.by,
				".RDS", sep="")
		}
		temp <- linker(no.taxa, arena.length, mean.log.individuals, length.parameter, 
			sd.parameter, max.distance, proportion.killed, competition.iterations, 
			no.plots, plot.length, concat.by, randomizations, cores,
			simulations, nulls, metrics)
		saveRDS(temp, file=filename)
	}
	return("Files saved to working directory")
}

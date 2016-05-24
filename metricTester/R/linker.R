#' Run spatial simulations, null and metric calculations to test metric + null performance
#'
#' This function wraps a number of wrapper functions into one big metric + null 
#' tester function. Only a single test is performed, with results saved into memory.
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
#' @param simulations Optional list of named spatial simulation functions to use. These
#' must be defined in the defineSimulations function. If invoked, this option will likely
#' be used to run a subset of the defined spatial simulations.
#' @param nulls Optional list of named null model functions to use. If invoked, this 
#' option will likely be used to run a subset of the defined null models.
#' @param metrics Optional list of named metric functions to use. If invoked, this option
#' will likely be used to run a subset of the defined metrics.
#' 
#' @details This function wraps a number of other wrapper functions into
#' one big metric + null performance tester function. Only a single test is performed, 
#' with results saved into memory. To perform multiple complete tests, use the
#' multiLinker function, which saves results to file.
#'
#' @return A list of lists of data frames. The first level of the output has one element 
#' for each simulation. The second level has one element for each null model. Each of
#' these elements is a list of two data frames, one that summarizes the plot-level
#' significance and another and arena-level significance.
#'
#' @export
#'
#' @importFrom geiger sim.bdtree
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' #below not run for timing issues on CRAN
#' #system.time(test <- linker(no.taxa=50, arena.length=300, mean.log.individuals=2, 
#' 	#length.parameter=5000, sd.parameter=50, max.distance=30, proportion.killed=0.2, 
#'	#competition.iterations=3, no.plots=15, plot.length=30, concat.by="richness", 
#'	#randomizations=3, cores="seq",
#'	#nulls=list("richness"=metricTester:::my_richnessNull,
#'	#"frequency"=metricTester:::my_frequency)))

linker <- function(no.taxa, arena.length, mean.log.individuals, length.parameter, 
	sd.parameter, max.distance, proportion.killed, competition.iterations, no.plots, 
	plot.length, concat.by, randomizations, cores, simulations, nulls, metrics)
{
	#set these things to NULL if they are not passed in, meaning that all defined sims,
	#nulls and metrics will be calculated
	if(missing(simulations))
	{
		simulations <- NULL
	}
	if(missing(nulls))
	{
		nulls <- NULL
	}
	if(missing(metrics))
	{
		metrics <- NULL
	}

	#simulate tree with birth-death process
	tree <- sim.bdtree(b=0.1, d=0, stop="taxa", n=no.taxa)
	#prep the data for spatial simulations
	prepped <- prepSimulations(tree, arena.length, mean.log.individuals, length.parameter, 
		sd.parameter, max.distance, proportion.killed, competition.iterations)
	#run the spatial simulations
	arenas <- runSimulations(prepped, simulations)
	#derive CDMs. plots are placed in the same places across all spatial simulations
	cdms <- multiCDM(arenas, no.plots, plot.length)
	
	#calculate observed metrics for all three spatial simulations
	observed <- lapply(cdms, function(x) observedMetrics(tree=tree, 
		picante.cdm=x$picante.cdm, metrics))
	#randomize all observed CDMs the desired number of times. this will generate a list of
	#lists of data frames. first level of list is for each spatial simulation (e.g. 3 if
	#there is random, habitat filtering and competitive exclusion). second level is for
	#randomizations, one for each. third level is data frames, one per null model
	allRandomizations <- lapply(1:length(cdms), function(x) metricsNnulls(tree=tree, 
		picante.cdm=cdms[[x]]$picante.cdm, regional.abundance=cdms[[x]]$regional.abundance,
		distances.among=cdms[[x]]$dists, cores=cores, 
		randomizations=randomizations, metrics=metrics, nulls=nulls))
	#reduce the randomizations to a list of lists of (first order of lists is for each
	#spatial simulation) data frames
	reduced <- lapply(allRandomizations, reduceRandomizations)
	#now lapply the errorChecker over each spatial simulation
	#the output of this is similar to above. list of lists of
	#data frames. first level for each simulation. second level for each null model.
	#the two data frames per second level summarizing the plot and arena-level
	#significance results
	results <- lapply(1:length(reduced), function(x) 
		errorChecker(observed=observed[[x]], reduced.randomizations=reduced[[x]],
		concat.by=concat.by, metrics))
	names(results) <- names(arenas)
	results
}

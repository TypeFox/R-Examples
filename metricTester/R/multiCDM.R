#' Wrapper for deriving CDMs from the results of multiple spatial simulations
#'
#' Given the results of a call to runSimulations(), this function places plots down
#' randomly (though identically across simulations). 
#'
#' @param simulations.result List of data frames of three columns: 
#' "individuals", "X", and "Y"
#' @param no.plots Number of plots to place
#' @param plot.length Length of one side of desired plot
#' 
#' @details Both the size and number of plots
#' are determined by the user. A conservative check (perhaps overly so) is in place to
#' ensure the function doesn't get stuck looking for solutions for how to randomly place
#' non-overlapping plots. Either decreasing the number or size of plots is
#' necessary if this throws an error.
#'
#' @return A list of data frames.
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
#' prepped <- prepSimulations(tree, arena.length=300, mean.log.individuals=2,
#' 	length.parameter=1000, sd.parameter=50, max.distance=20, proportion.killed=0.2,
#'	competition.iterations=5)
#'
#' #run the spatial simulations
#' arenas <- runSimulations(prepped)
#'
#' #derive CDMs. plots are placed in the same places across all spatial simulations.
#' #density of individuals per arena is low enough in this example that sometimes all
#' #plots contain < 2 species, and are cut, causing an error. this not run so as not to
#' #throw errors on CRAN
#' #cdms <- multiCDM(arenas, no.plots=10, plot.length=20)

multiCDM <- function(simulations.result, no.plots, plot.length)
{
	results <- lapply(1:length(simulations.result), function(x)
		makeCDM(single.simulation=simulations.result[[x]], no.plots=no.plots,
		plot.length=plot.length))
	names(results) <- names(simulations.result)
	results
}

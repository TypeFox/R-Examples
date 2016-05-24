#' Prep data for spatial simulations
#'
#' Given the required parameters for defined spatial simulations, will prepare an
#' object of class simulations.input for actual simulation.
#'
#' @param tree Phylo object
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
#' 
#' @details This function preps the input for any of the spatial simulations as defined in
#' defineSimulations. If additional parameters are ever required for those simulations,
#' they would have to be added as additional arguments here.
#'
#' @return A prepared simulations.input object
#'
#' @export
#'
#' @importFrom stats pnorm quantile rlnorm rnorm runif
#' 
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' tree <- geiger::sim.bdtree(b=0.1, d=0, stop="taxa", n=50)
#'
#' prepped <- prepSimulations(tree, arena.length=300, mean.log.individuals=2, 
#' 	length.parameter=5000, sd.parameter=50, max.distance=20, proportion.killed=0.2,
#'	competition.iterations=3)

prepSimulations <- function(tree, arena.length, mean.log.individuals, length.parameter, 
	sd.parameter, max.distance, proportion.killed, competition.iterations)
{
	dat <- list("tree"=tree, "arena.length"=arena.length, 
	"mean.log.individuals"=mean.log.individuals, "length.parameter"=length.parameter, 
	"sd.parameter"=sd.parameter, "max.distance"=max.distance, 
	"proportion.killed"=proportion.killed, "competition.iterations"=competition.iterations)
	class(dat) <- c("list", "simulations.input")
	dat
}

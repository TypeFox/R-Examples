#' Plot simulated plots in arena
#'
#' Given a matrix of plot bounds, plots the plots in an already plotted,
#' simulated arena.
#'
#' @param plot.bounds Matrix of plot bounds
#' 
#' @details Plots plots as defined by the supplied matrix, e.g. a call to
#' plotPlacer. An active plot with the simulated arena needs to already be open, 
#' see example.
#'
#' @return Plotted plots
#'
#' @export
#'
#' @importFrom graphics plot polygon
#' @importFrom colorRamps blue2green2red
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#'
#' tree <- geiger::sim.bdtree(b=0.1, d=0, stop="taxa", n=50)
#'
#' temp <- evolveTraits(tree)
#'
#' phydistmatrix <- ape::cophenetic.phylo(temp[[1]])
#'
#' #define a color for each species
#' cols <- colorRamps::blue2green2red(nrow(phydistmatrix))
#'
#' #prep the data for the simulation
#' prepped <- prepSimulations(tree, arena.length=300, mean.log.individuals=2, 
#' length.parameter=5000, sd.parameter=50, max.distance=20, proportion.killed=0.2,
#' competition.iterations=3)
#'
#' positions <- filteringArena(prepped)
#'
#' #plot the arena. don't close the window
#' plot(positions$arena$X, positions$arena$Y, pch=20, cex=0.5, xlim=c(0,300), ylim=c(0,300), 
#' col=cols[positions$arena$individuals])
#'
#' bounds <- plotPlacer(no.plots=10, arena.length=300,
#'	plot.length=50)$plot.bounds
#'
#' plotPlotter(bounds)

plotPlotter <- function(plot.bounds)
{
	for(i in 1:dim(plot.bounds)[1])
	{
		polygon(c(plot.bounds[i,1],plot.bounds[i,2],plot.bounds[i,2],
		plot.bounds[i,1]),c(plot.bounds[i,3],plot.bounds[i,3],
		plot.bounds[i,4],plot.bounds[i,4]))
	}
}

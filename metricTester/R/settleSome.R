#' Randomly settle individuals in a spatial arena
#'
#' Given output from the killSome function, randomly settles individuals in the arena.
#'
#' @param killSome.output Output from the killSome function
#' 
#' @details This function uses the number killed element of the killSome output to
#' randomly draw from the regional abundance vector, then settles the individuals at
#' random in the arena. 
#'
#' @return A list of 4 elements: the average relatedness in the geographic neighbordhood
#' of consideration (passed directly from the killSome output, not re-calculated here), 
#' the regional abundance vector, the new spatial arena, and the dimensions of that arena.
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
#' #prep the data for the simulation
#' prepped <- prepSimulations(tree, arena.length=300, mean.log.individuals=2, 
#' length.parameter=5000, sd.parameter=50, max.distance=20, proportion.killed=0.2,
#' competition.iterations=5)
#'
#' #use the competition simulation
#' positions <- competitionArena(prepped)
#'
#' #in normal use, these parameters will be carried down from the simulations.input object
#' new.arena <- killSome(tree, arena.output=positions, max.distance=50, 
#' proportion.killed=0.2)
#'
#' #now settle some indiviudals
#' newer.arena <- settleSome(new.arena)
#'
#' #look at how number of individuals in arena changes
#' dim(new.arena$arena)
#' dim(newer.arena$arena)

settleSome <- function(killSome.output)
{
	#sample the same number of individuals you killed from the regional abundance vector
	individuals <- sample(killSome.output$regional.abundance, 
		size=killSome.output$no.killed)
	
	#start a dataframe to bind X,Y coordinates into	
	to.bind <- data.frame(individuals)

	#generate random X,Y coordinates that are evenly distributed across the arena
	to.bind$X <- runif(n=length(individuals), min=0, max=killSome.output$dims[2])
	to.bind$Y <- runif(n=length(individuals), min=0, max=killSome.output$dims[2])
	
	output <- list(related=killSome.output$related, 
		regional.abundance=killSome.output$regional.abundance, 
		arena=rbind(killSome.output$arena, to.bind), dims=killSome.output$dims)
}

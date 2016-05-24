#' Simulate competitive exclusion over multiple generations
#'
#' Given a simulations.input object, creates an initial random arena, then removes some
#' closely related individuals in the arena, settles individuals based on abundances from
#' a regional species pool, and repeats across the desired number of generations.
#'
#' @param simulations.input A prepared simulations.input object from prepSimulations
#' 
#' @details This function combines the killSome and settleSome functions into a loop that
#' runs for the desired number of generations. If the number of individuals in the arena
#' exceeds 50,000, then the killSomeBig function is invoked. This runs much slower than
#' the default version of the function, but will not use all memory and crash at large
#' numbers of individuals.
#'
#' @return A list of 5 elements: the average relatedness in the geographic neighborhood
#' of consideration (appended to any previous values that were fed into the function), 
#' the number of individuals killed, the original input regional
#' abundance vector, the new spatial arena, and the dimensions of that arena. On the last
#' iteration, it returns the arena BEFORE settling new individuals randomly.
#'
#' @export
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
#' 	competition.iterations=2)
#'
#' test <- competitionArena(prepped)

competitionArena <- function(simulations.input)
{
	initialArena <- randomArena(simulations.input)

	for(i in 1:simulations.input$competition.iterations)
	{
		#take the initialArena and kill off some of the individuals in genetically 
		#clustered neighborhoods. the cutoff to invoke killSomeBig should be somewhere
		#around 50,000 individuals
		if(dim(initialArena$arena)[1] < 50000)
		{
			killed.arena <- killSome(simulations.input$tree, initialArena, 
				simulations.input$max.distance, simulations.input$proportion.killed)
		}
		else
		{
			killed.arena <- killSomeBig(simulations.input$tree, initialArena, 
				simulations.input$max.distance, simulations.input$proportion.killed)
		}

		#add individuals back in, but save this as initialArena, so that it gets plugged
		#back in next iteration instead of the original random arena
		initialArena <- settleSome(killed.arena)
	}

	return(killed.arena)
}

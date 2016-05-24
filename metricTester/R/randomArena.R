#' Generate a random spatial arena
#'
#' Given a simulations.input object, will create a randomly settled arena according to the
#' parameters in the simulations.input object.
#'
#' @param simulations.input A prepared simulations.input object from prepSimulations
#' 
#' @details This function uses the log-normal regional abundance distribution to randomly
#' assign abundances to species in the tree. It then draws from this regional abundance
#' distribution to settle individuals at random in the landscape. 
#'
#' @return A list of 4 elements: the mean of the genetic distance matrix of the input
#' phylogeny, the regional abundance vector (where each element is a species name, 
#' repeated as many times as is present in pool), the spatial arena, and the dimensions of
#' that arena.
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
#' length.parameter=5000, sd.parameter=50, max.distance=20, proportion.killed=0.2,
#' competition.iterations=2)
#'
#' test <- randomArena(prepped)

randomArena <- function(simulations.input)
{
	#this ugly bit of code is because I used to define these things as arguments to funct
	#this makes input easier
	x.min=0
	x.max=simulations.input$arena.length
	y.min=0
	y.max=simulations.input$arena.length

	#generate log-normal regional abundance curve, and randomly assign abundances to species
	indivs.per.species <- rlnorm(n=length(simulations.input$tree$tip.label), 
		simulations.input$mean.log.individuals, sdlog=1)
	
	#set species with < 0 individuals to 0 abundance
	indivs.per.species[indivs.per.species < 0] <- 0

	#round abundances to no decimal places
	indivs.per.species <- round(indivs.per.species)

	#actually generate a vector individuals with species identities (the "regional pool")
	individuals <- c()

	individuals <- rep(simulations.input$tree$tip.label, times=indivs.per.species)

	#start a dataframe to bind X,Y coordinates into	
	arena <- data.frame(individuals)

	#generate random X,Y coordinates centered around the middle of the arena
	arena$X <- sample(x.min:x.max, size=length(individuals), replace=TRUE)
	arena$Y <- sample(y.min:y.max, size=length(individuals), replace=TRUE)
	
	#create and return the output
	
	output <- list(related=mean(ape::cophenetic.phylo(simulations.input$tree)),
		regional.abundance=individuals, arena=arena, 
		dims=c(x.min, x.max, y.min, y.max))
	
	return(output)
}

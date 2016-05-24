#' Simulate a community assembled according to habitat filtering
#'
#' Given a simulations.input object, will create an arena settled according to habitat
#' filtering rules (and parameters defined by prepSimulations).
#'
#' @param simulations.input A prepared simulations.input object from prepSimulations
#' 
#' @details This simulation has been updated to avoid a previous consequence where
#' individuals clumped near the center of the arena. Now, species' spatial preferences,
#' which previously would have approximated a normal distribution around the center of the
#' arena, are smoothed to a uniform distribution. Thus, individuals are more evenly
#' distributed throughout the simulated arena.
#'
#' @return A list of 3 elements: the original input regional
#' abundance vector, the new spatial arena, and the dimensions of that arena. 
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
#' test <- filteringArena(prepped)

filteringArena <- function(simulations.input)
{
	#evolve two traits independently up phylogeny. returns a list where first element is
	#input tree and second element is a matrix of traits (x,y spatial preferences)
	evolveTraits.results <- evolveTraits(simulations.input$tree)

	#IN THE ORIGINAL SIMULATIONS, WE DID NOT RUN EITHER OF THESE NEXT TWO LINES OF DATA
	#TO TRANSFORM TO UNIFORM DISTRIBUTION
	
	#determine what quantile first trait falls in and reclassify the trait to that
	#this essentially turns the normal distribution uniform
	evolveTraits.results[[2]][,1] <- pnorm(evolveTraits.results[[2]][,1], 
		mean(evolveTraits.results[[2]][,1]), sd(evolveTraits.results[[2]][,1]))
	
	#do same for second trait
	evolveTraits.results[[2]][,2] <- pnorm(evolveTraits.results[[2]][,2], 
		mean(evolveTraits.results[[2]][,2]), sd(evolveTraits.results[[2]][,2]))

	#scale results to size of arena
	scaled.results <- scaler(evolveTraits.results[[2]], min.arena=0, 
		max.arena=simulations.input$arena.length)

	#simulate the number of individuals per species to follow a log-normal abundance dist
	indivs.per.species <- rlnorm(n=length(evolveTraits.results[[1]]$tip.label),
		simulations.input$mean.log.individuals, sdlog=1)
	indivs.per.species[indivs.per.species < 0] <- 0

	indivs.per.species <- round(indivs.per.species)

	individuals <- c()

	#expand out the abundance distributions to actual identities
	for(i in 1:length(indivs.per.species))
	{
		individuals <- append(individuals, rep(evolveTraits.results[[1]]$tip.label[i], 
			times=indivs.per.species[i]))
	}
	
	arena <- data.frame(individuals)

	X <- c()
	Y <- c()
	
	#per species simulate a normal distribution with the length and sd parameters provided
	#around each species spatial preference. 
	for(i in 1:length(individuals))
	{
		X.options <- rnorm(n=simulations.input$length.parameter, 
			mean=scaled.results[row.names(scaled.results)==individuals[i], 1], 
			sd=simulations.input$sd.parameter)
		X[i] <- sample(X.options, size=1)

		Y.options <- rnorm(n=simulations.input$length.parameter, 
			mean=scaled.results[row.names(scaled.results)==individuals[i], 2], 
			sd=simulations.input$sd.parameter)
		Y[i] <- sample(Y.options, size=1)
	}
	
	arena$X <- X
	arena$Y <- Y

	#this ugly piece of code for compatibility with older pieces of code I should revise	
	x.min=0
	x.max=simulations.input$arena.length
	y.min=0
	y.max=simulations.input$arena.length
	dims=c(x.min, x.max, y.min, y.max)
	
	results <- list("regional.abundance"=as.character(arena$individuals), "arena"=arena, 
		"dims"=dims)
	
	return(results)
}

#' Remove most closely related individuals
#'
#' Given a phylogenetic tree, a spatial arena of individuals with species identities,
#' and arguments for the desired distance and percent removed, removes some of the most
#' closely related individuals in the arena.
#'
#' @param tree Phylo object
#' @param arena.output A spatial arena with three columns: individuals (the species ID), 
#' X (the x axis location of that individual), and Y (the y axis location). The
#' arena.output actually needs a number of other elements in order for later functions to
#' work properly, so any modifications to the code should take note of this.
#' @param max.distance The geographic distance within which geographically neighboring
#' indivduals should be considered to influence the individual in question.
#' @param proportion.killed The percent of individuals in the total arena that should be
#' considered (as a proportion, e.g. 0.5 = half).
#' 
#' @details This function identifies individuals in the most genetically clustered
#' geographic neighborhoods, continues on to identify the most closely related individual
#' to a focal individual, and randomly chooses whether to remove that individual or the
#' focal individual. It expects a list with a number of additional elements beyond the 
#' arena (currently, the mean genetic relatedness of geographic neighborhoods, a vector of
#' regional abundance [where each element is a species name, repeated as many times as is
#' present in pool], and the dimensions of the arena). 
#'
#' @return A list of 5 elements: the average relatedness in the geographic neighbordhood
#' of consideration (appended to any previous values that were fed into the function), 
#' the number of individuals killed, the original input regional
#' abundance vector, the new spatial arena, and the dimensions of that arena.
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
#' 	length.parameter=5000, sd.parameter=50, max.distance=20, proportion.killed=0.2,
#' 	competition.iterations=5)
#'
#' #use the competition simulation
#' positions <- competitionArena(prepped)
#'
#' #in normal use, these parameters will be carried down from the simulations.input object
#' new.arena <- killSome(tree, arena.output=positions, max.distance=50, 
#' 	proportion.killed=0.2)
#'
#' #look at how number of individuals in arena changes
#' dim(positions$arena)
#' dim(new.arena$arena)

killSome <- function(tree, arena.output, max.distance, proportion.killed)
{
	#set up a blank output list to save into
	output <- list()

	#save the species identities of all individuals
	individual.identities <- arena.output$arena$individuals

	#create a genetic distance matrix
	gen.dist <- ape::cophenetic.phylo(tree)
	
	#create a matrix of individuals for use in geographic distance calculations. obviously 
	#very similar to the input data frame, but dist doesn't work right with data frames
	for.geo.dist <- matrix(cbind(arena.output$arena$X, arena.output$arena$Y), ncol=2)
	
	#calculate all pairwise geographic distances
	geo.dist <- dist(for.geo.dist)
	
	#convert the distance matrix into a real matrix
	geo.dist.matrix <- as.matrix(geo.dist)

	#set all geographic distances from one individual to itself to NA	
	diag(geo.dist.matrix) <- NA
	
	#set all geographic distances greater than the max distance to consider to NA
	geo.dist.matrix[geo.dist.matrix > max.distance] <- NA
	
	#replace these geographic distances with genetic distances, so we can calculate the 
	#genetic neighborhood of individuals. first figure out who the species involved in
	#each pairwise comparison are
	
	#the individuals involved in each comparison
	individual.involved <- which(!is.na(geo.dist.matrix), arr.ind=TRUE)
	
	#replace with species' identities. there is probably a better way of doing this, but
	#subsetting in this manner seems to work, and just spins right through the whole data
	#frame one column of individual.involved at a time. if i split the results in two and
	#stack into two columns it gives the right results
	species.involved <- matrix(individual.identities[individual.involved], 
		nrow=length(individual.identities[individual.involved])/2, ncol=2)
	
	#find the appropriate genetic distances for all these species combinations
	specific.gen.dist <- gen.dist[species.involved]
	
	#plug these genetic distances into the geographic distance matrix of just the closest 
	#individuals
	specific.gen.dist.matrix <- geo.dist.matrix
	specific.gen.dist.matrix[which(!is.na(specific.gen.dist.matrix))] <- specific.gen.dist
	
	#find the average relatedness of every individual to its closest neighbors
	average.relatedness <- apply(specific.gen.dist.matrix, 2, mean, na.rm=TRUE)
	
	###DELETE THIS CODE AFTER TESTING. you just want to know mean average relatedness in 
	#geographic neighborhoods with iterations
	
	output$related <- append(arena.output$related, mean(average.relatedness))
	
	#define the mean genetic distance below which individuals will be considered to be in 
	#genetically clustered geographic neighborhoods
	cutoff <- quantile(average.relatedness, probs=proportion.killed, na.rm=TRUE)
	
	#find the individual IDs of those in the most genetically clustered neighborhoods
	individuals.considered <- which(average.relatedness <= cutoff)
	
	####NEED TO INCLUDE AN ARGUMENT HERE TO NOT CONSIDER MIN of ZERO IF YOU WANT TO NOT
	####CONSIDER CONSPECIFICS

	#find the minimum genetic distance between every individual (that has a geographic
	#neighbor within the maximum distance), and its nearest genetic neighbor
	mins <- suppressWarnings(apply(specific.gen.dist.matrix, 2, min, na.rm=TRUE))
	
	#make a temporary table where you bind the minimum value onto the bottom of the table
	temp.table <- rbind(specific.gen.dist.matrix, mins)
	
	#use the compareMins utility function to compare the last element in a column to all
	#other elements in that column
	semifinal.matrix <- apply(temp.table, 2, compareMins)
	
	final.matrix <- semifinal.matrix
	
	#set all FALSE elements in the table to NA
	final.matrix[final.matrix==FALSE] <- NA
	
	#then identify who the most closely related individual is to every other individual
	closest.table <- which(!is.na(final.matrix), arr.ind=TRUE)
	
	#subset this table to only those individuals that were in genetically clustered 
	#geographic neighborhoods
	final.table <- closest.table[closest.table[,2] %in% individuals.considered,]
	
	#randomly select half of the individuals to kill. BE SURE TO ONLY TAKE THE UNIQUE
	#individuals involved in each comparison, or your plots will get denser with time
	kill.list <- unique(sample(final.table[,1], size=(dim(final.table)[1])/2))

	new.arena <- arena.output$arena[-kill.list,]
	
	#create and return the output
	
	output$no.killed=length(kill.list)
	output$regional.abundance=arena.output$regional.abundance
	output$arena=new.arena
	output$dims=arena.output$dims
	
	return(output)
}

#' Remove most closely related individuals for large arenas
#'
#' Given a phylogenetic tree, a large spatial arena of individuals with species
#' identities, and arguments for the desired distance and percent removed, removes some of
#' the most closely related individuals in the arena.
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
#' present in pool], and the dimensions of the arena). This function is invoked if the
#' number of individuals in the arena exceeds 50,000. It was originally intended only for 
#' arenas with large numbers of individuals, but it ends up being faster than the original
#' version of the function for all but the smallest arenas. That said, if the local area
#' within the arena that the function considers does not contain any individuals, then it
#' will fail. To overcome this, either a larger max.distance would need to be set, or more
#' individuals should be settled in the arena, or the original version of the function
#' should be invoked by modifying the code in competitionArena(). Note, however, that if 
#' there are more than e.g. 50,000 individuals total in the arena, the original version of
#' the function should never be used, or all memory in the computer will be used and it
#' will crash.
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
#' prepped <- prepSimulations(tree, arena.length=300, mean.log.individuals=1, 
#' length.parameter=5000, sd.parameter=50, max.distance=20, proportion.killed=0.2,
#' competition.iterations=5)
#'
#' #use the random simulation to generate a starting arena
#' random <- randomArena(prepped)
#'
#' #in normal use, these parameters will be carried down from the simulations.input object
#' new.arena <- killSomeBig(tree, arena.output=random, max.distance=50, 
#' proportion.killed=0.2)
#'
#' #look at how number of individuals in arena changes
#' dim(random$arena)
#' dim(new.arena$arena)

killSomeBig <- function(tree, arena.output, max.distance, proportion.killed)
{
	#set up a blank output list to save into
	output <- list()

	#create a genetic distance matrix
	genDists <- ape::cophenetic.phylo(tree)

	#repeatedly sample local neighborhoods to determine what distribution of genetic
	#distances among geographic neighbors is. first figure out how many plots you will
	#place each iteration of this scheme. the math from this comes from the guestimate in
	#plotPlacer(). basically, solve for what the number of plots would be if you
	#wanted to sample an area that was 40% of the total arena, then round down
	
	noPlots <- 0.4 * ((max(arena.output$dims)^2)/(max.distance^2))
	noPlots <- floor(noPlots)
	
	avRelatedness <- list()
	
	#arbitrarily choose 10 loops of sampling and calculating MPD	
	for(i in 1:10)
	{
		temp <- makeCDM(single.simulation=arena.output, no.plots=noPlots,
			plot.length=max.distance)
		avRelatedness[[i]] <- modifiedMPD(samp=temp$picante.cdm, dis=genDists, 
			abundance.weighted="interspecific")
	}

	#determine what the cutoff is below which approximately proportion.killed of
	#individuals should be considered as in genetically clustered geographic neighborhoods
	avRelatedness <- unlist(avRelatedness)
	cutoff <- quantile(avRelatedness, proportion.killed)

	#go into a for loop here, since cannot calculate a geographic distance matrix among
	#large numbers of individuals. the general plan to circumvent that is to calculate
	#the distance among a smaller subset of individuals as defined by max.distance. if
	#there are still too many individuals in these smaller segments, return an error
	
	#first make a temporary data frame where individuals get unique identifying numbers
	#and randomize the order of this so it does not proceed from s1 to 2, etc.
	tempArena <- arena.output$arena
	tempArena$identifier <- 1:(dim(arena.output$arena)[1])
	tempArena <- tempArena[sample(row.names(tempArena)), ]

	for(i in 1:dim(tempArena)[1])
	{
		#subset the arena to just those individuals within max.distance (in x,y sense) of
		#focal individual. first get its X, Y location. if you try with speed.factor, e.g.
		#focalY <- tempArena$Y[tempArena$identifier %in% individualsConsidered[i]]
		focalX <- tempArena$X[i]
		focalY <- tempArena$Y[i]
		
		#it is possible that the focal individual could have been killed on a previous
		#iteration. add a check to skip to the next iteration if so
		if(length(focalX) < 1)
		{
			next()
		}
		
		else
		{
			minX <- focalX - max.distance
			maxX <- focalX + max.distance
			minY <- focalY - max.distance
			maxY <- focalY + max.distance
			miniArena <- tempArena[tempArena$X >= minX,]
			miniArena <- miniArena[miniArena$X <= maxX,]
			miniArena <- miniArena[miniArena$Y >= minY,]
			miniArena <- miniArena[miniArena$Y <= maxY,]
		}
		
		#add a check here so that if this miniArena still has too many individuals in it
		#it will throw an error before running out of memory trying to calculate distances
		if(dim(miniArena)[1] > 10000)
		{
			stop("You are considering too many individuals and would run out of memory")
		}
		#note that you might want to add an else if statement here that if miniArena does
		#not have any individuals in it, stop and modify some arguments. otherwise this
		#function will fail (the specific error it throws comes with the min(focalGenDist)
		#below). i do not currently have any problems with this, but killSomeBig gets
		#invoked on an absolute number of individuals level. it will fail on a density of
		#individuals in the arena level though. so if you make a really big but sparse
		#arena the whole thing will fail
		else
		{
			#calculate the pairwise distances among all individuals and format as a matrix
			geoDists <- as.matrix(dist(miniArena[,2:3]))
			#give row and col names to these to ensure you tabulate things properly 
			#(I do not think this is actually necessary)
			row.names(geoDists) <- miniArena$identifier
			colnames(geoDists) <- miniArena$identifier
		}
		
		#subset to just the distances between the focal individual and others
		specDists <- geoDists[rownames(geoDists) == tempArena$identifier[i],]

		#exclude any distances that are larger than max.distance (remember we only did it
		#in an xy sense before)
		specDists <- specDists[specDists <= max.distance]

		#subset miniArena to these geographically close individuals
		toSample <- miniArena[miniArena$identifier %in% names(specDists),]

		#calculate the interspecific MPD of these individuals and retain to decide which
		#individuals to actually engage in competition. only consider individuals in
		#neighborhoods more clustered than the cutoff defined above. note that individuals
		#can be NA (e.g. only individual in that neighborhood), so need to add is.na below
		indivRelated <- modifiedMPD(samp=t(summary(toSample$individuals)),
			dis=genDists)
		
		if(indivRelated > cutoff | is.na(indivRelated))
		{
			next()
		}
		
		else
		{
			#exclude the geographic distance to the individual in question and proceed
			specDists <- specDists[names(specDists) != tempArena$identifier[i]]
		}

		#now determine the species identities of the individuals being considered
		individualIdentities <- miniArena$individuals[miniArena$identifier 
			%in% names(specDists)]
		
		#subset the entire genetic distance matrix to the row for the focal species
		tempGenDists <- genDists[row.names(genDists) == tempArena$individuals[i],]
	
		#subset this to genetic distances between the focal individual and all possible
		#species from individualIdentities. NOTE THAT YOU DO NOT WANT TO EXCLUDE THE 
		#genetic distance of 0, as it is possible that another individual is same species
		#and you want to consider these. note, it is also possible that the focal indiv
		#is equally closely related to different species, so sample from most closely
		#related species
		focalGenDists <- tempGenDists[names(tempGenDists) %in% individualIdentities]
		
		minGenDist <- min(focalGenDists)
		
		#find the species identity of the most closely related species in the geographic
		#neighborhood under consideration (randomly sampling among closely related spp)
		#you will kill off an individual of this species
		
		#first, need to account for the fact that branches are ever so slightly different
		#lengths, even though the tree is ultrametric. this returns the names of any
		#species that are equally closely related to the species that ape/R thinks is the
		#most closely related
		
		possibleSp <- which(sapply(focalGenDists, function(x) 
			all.equal(minGenDist, x, check.names=FALSE))==TRUE)
		
		toKill <- sample(names(possibleSp), 1)
		
		#determine identifier of an individual of the species toKill within the 
		#max.distance of the individual considered, and which is not the focal individual
		toSample <- toSample[toSample$identifier != tempArena$identifier[i],]

		#just "randomly" take the first row for simplicity sake
		sampled <- toSample[toSample$individuals == toKill,][1,]
		
		#remove this individual from the output.arena (and from successively smaller
		#versions of the arena as you remove individuals)
		tempArena <- tempArena[tempArena$identifier != sampled$identifier,]
	}

	#get rid of the "identifier" column from tempArena. note that at the moment, this 
	#function kills fewer individuals than you input. that is because certain
	#closely related individuals get killed off at early stages of the for loop, then
	#they are already taken out of tempArena and the last line of the for loop does
	#nothing
	
	tempArena$identifier <- NULL

	#make a CDM so that you can calculate local relatedness within the arena. 
	#note that makeCDM needs a list with elements arena and dims, so make a dummy list
	
	tempArena2 <- list()
	tempArena2$arena <- tempArena
	tempArena2$dims <- arena.output$dims
	
	tempCDM <- makeCDM(single.simulation=tempArena2, no.plots=noPlots, 
		plot.length=max.distance)
	
	#calculate the interspecific MPD for these plots
	tempRelated <- modifiedMPD(samp=tempCDM$picante.cdm, dis=genDists, 
		abundance.weighted="interspecific")
	
	#calculate the mean interspecific MPD for the plots, and append to vector in output
	related <- mean(tempRelated)
	
	output$related <- append(arena.output$related, related)

	output$no.killed=length(arena.output$regional.abundance) - length(tempArena$X)
	output$regional.abundance=arena.output$regional.abundance
	output$arena=tempArena
	output$dims=arena.output$dims
	
	return(output)
}

#' Randomize community data matrix with dispersal null model
#'
#' Null model that maintains species richness and settles species in simulated plots 
#' with a probability proportional to their abundance in and distance from the plot in
#' question.
#'
#' @param picante.cdm Picante-style community data matrix with
#' communities/plots/plots/etc as rows and species as columns
#' @param tree Phylo object
#' @param distances.among A symmetric distance matrix, summarizing the distances among all
#' plots from the cdm.
#' @param abundance.matters Default is TRUE. If FALSE, species are sampled from
#' neighboring grid cells with equal probability.
#' @param abundance.assigned Default is "directly". If set to "explore", rather than
#' assigning species abundances equal to those in neighboring cells, a normal distribution
#' of standard deviation 1 is created around that abundance, rounded to whole numbers,
#' negative numbers are converted to 1, and then abundance is assigned from this
#' distribution. Note, importantly, that if the cdm is a relative abundance matrix, as
#' opposed to an absolute abundance matrix, if "explore" is chosen, odd results are likely
#' (don't do it). If set to "overall", which may be preferable when abundance.matters is
#' set to FALSE, species' are assigned abundances by drawing from the vector of non-zero
#' abundances from the original cdm.
#' 
#' @details This function can run quite slowly, as it employs a while loop to discard
#' selected species if they are already contained in the plot being simulated. The
#' function contains an internal check to ensure that it doesn't get stuck in an
#' indefinite while loop. Specifically, it checks that an observed plot does not
#' contain more species than the sum of unique species in the remaining plots. The
#' argument distances.among is flexible, and can relate to e.g., straight-line distances,
#' great-circle distances, and ecological distances. If the argument abundance.matters is
#' set to FALSE, then species' abundances in neighboring grid cells does not influence
#' their probability of settling (only their proximity influences the probability of
#' settling). If this is the case, it is recommended that the argument abundance.assigned
#' be set to "overall".
#'
#' @return A matrix with the same dimensions as the input cdm.
#'
#' @export
#'
#' @importFrom picante sample2matrix
#'
#' @references Miller, E. T. 2016. A new dispersal-informed null model for
#' community ecology shows strong performance. biorxiv.
#'
#' @examples
#' #set up a matrix to simulate lat/long
#' coordDF <- matrix(ncol=2, nrow=100)
#'
#' coordDF[,1] <- runif(n=100, min=40, max=50)
#' coordDF[,2] <- runif(n=100, min=-130, max=-120)
#'
#' #convert to data frame, give column names. also give row names such as if the cells had
#' #names (as they should or there'd be no way to track them)
#' coordDF <- as.data.frame(coordDF)
#'
#' row.names(coordDF) <- paste("cell", 1:100, sep="")
#'
#' names(coordDF) <- c("latitude","longitude")
#'
#' #calculate the distances among all of these points. in the real program you're going to
#' #want to calculate great arc distance or whatever it's called
#' distances <- dist(coordDF, diag=TRUE, upper=TRUE)
#'
#' #turn it into a symmetric distance matrix
#' distances <- as.matrix(distances)
#'
#' #simulate a regional phylogeny of 100 species
#' tree <- geiger::sim.bdtree(b=1, d=0, stop="taxa", n=100)
#'
#' #simulate a community data matrix of 100 cells by 100 species. do it 4 times so that
#' #you can use your simulateComm function and have it span a reasonable range of richness
#' sim.abundances <- round(rlnorm(5000, meanlog=2, sdlog=1)) + 1
#'
#' cdm1 <- simulateComm(tree, richness.vector=10:34, abundances=sim.abundances)
#' cdm2 <- simulateComm(tree, richness.vector=10:34, abundances=sim.abundances)
#' cdm3 <- simulateComm(tree, richness.vector=10:34, abundances=sim.abundances)
#' cdm4 <- simulateComm(tree, richness.vector=10:34, abundances=sim.abundances)
#'
#' #bind these into a list and use dplyr rbind_all to bind together. recast as data frame
#'
#' cdmList <- list(cdm1, cdm2, cdm3, cdm4)
#'
#' cdm <- dplyr::rbind_all(cdmList)
#'
#' cdm <- as.data.frame(cdm)
#'
#' #fix as necessary manually here (i.e. make sure dimensions are 100 x 100), seems to 
#' #usually work. then give cell names
#'
#' row.names(cdm) <- paste("cell", 1:100, sep="")
#'
#' #fill NAs with 0s.
#'
#' cdm[is.na(cdm)] <- 0
#'
#' #not run
#' #newCDM <- dispersalNull(cdm, tree, distances)

dispersalNull <- function(picante.cdm, tree, distances.among, abundance.matters=TRUE, 
	abundance.assigned="directly")
{
	#create a check to ensure that no cell in the CDM contains more species than are found
	#in the sum of the remaining cells
	tempCheck <- checkCDM(picante.cdm)
	if(tempCheck == "fail")
	{
		stop("CDM incompatible with dispersalNull model. See 'checkCDM' for details")
	}

	#create a check to ensure that all species that occur in the cdm are also in the tree
	if(length(setdiff(names(picante.cdm),tree$tip.label))!=0)
	{
		stop("You have included species in your cdm that are not in your phylogeny")
	}
	
	#create a check to ensure that cdm and distances.among contain the same exact cells
	if(length(setdiff(row.names(picante.cdm), row.names(distances.among))) != 0 &
		length(setdiff(row.names(distances.among), row.names(picante.cdm))) != 0)
	{
		stop("Your cdm plot names and distance matrix names do no match")
	}
	
	#create a check to ensure that cdm and distances.among are in the same order cell-wise
	if(any(row.names(picante.cdm)!=row.names(distances.among)))
	{
		stop("Your cdm and distance matrix are not in the same plot order")
	}

	#calculate species richness of input cdm. for each plot, sample a grid cell based
	#on its proximity to the plot in question (this is your sampleNear function), then 
	#sample a species based on its abundance in the sampled grid cell. if that species is
	#already included in the replacement plot, run this over again. do this entire loop
	#per plot while the length of replacement species is < the observed richness for
	#that plot
	
	richness <- apply(picante.cdm, 1, lengthNonZeros)

	#generate a vector of all non-zero abundances from the input CDM
	overallAbundance <- picante.cdm[picante.cdm!=0]

	#create a list that you will use to save vectors of species into (one list element
	#per plot)
	replacementList <- list()
	
	for(i in 1:dim(picante.cdm)[1])
	{
		#create a temporary empty data frame in phylocom format. make it just long enough
		#for the plot in question
		phylocom <- matrix(ncol=3, nrow=richness[i], 0)
		phylocom <- as.data.frame(phylocom)
		#give it appropriate names
		names(phylocom) <- c("plot", "abund", "id")
		#start the while loop. set j equal to 0. this will be like a row ID. each time
		#an appropriate species is found, it will add one to the ID and move along the
		#phylocom df. the while loop continues while the length of species for a given
		#plot remains less than the richness of the observed plot
		j <- 0
		while(length(phylocom[phylocom$plot==row.names(picante.cdm)[i],]$id) < richness[i])
		{
			#select the plot you will sample from. give your selectNear a column from
			#the distance matrix
			selectedPlot <- selectNear(distances.among[,i])
			#select a species from that plot. probability proportional to abundance
			if(abundance.matters)
			{
				temp <- sample(x=picante.cdm[selectedPlot,], size=1, 
					prob=picante.cdm[selectedPlot,])
			}
			#or where probability is not proportional to abundance in neighboring cell.
			#note that you need to restricted the species from the selected plot to be
			#those with non-zero abundances! previously did not, this was a bug
			else
			{
				#note that this means the selected species lose their names, so add them
				#back in
				possible <- picante.cdm[selectedPlot,][picante.cdm[selectedPlot,] !=0 ]
				names(possible) <- names(picante.cdm[selectedPlot,])[picante.cdm[selectedPlot,]!=0]
				temp <- sample(x=possible, size=1)
			}
			#if the species selected is not found in that plot in the growing phylocom
			#data frame, add the relevant info to that frame
			if(!(names(temp) %in% phylocom[phylocom$plot==row.names(picante.cdm)[i],]$id))
			{
				#move the row ID along
				j <- j+1
				#set the plot name to be correct
				phylocom[j,1] <- row.names(picante.cdm)[i]
				#assign the species the abundance it had in the original CDM
				if(abundance.assigned=="directly")
				{
					phylocom[j,2] <- temp
				}
				#or else create a normal distribution and sample from it
				else if(abundance.assigned=="explore")
				{
					distribution <- round(rnorm(n=100, mean=as.numeric(temp), sd=1))
					distribution[distribution < 0] <- 1
					chosen <- sample(distribution, 1)
					phylocom[j,2] <- chosen
				}
				else if(abundance.assigned=="overall")
				{
					chosen <- sample(overallAbundance, 1)
					phylocom[j,2] <- chosen
				}
				else
				{
					stop("abundance.assigned argument set to unrecognized value")
				}
				phylocom[j,3] <- names(temp)
			}
		}
		#set the appropriate element in the list to the temporary dataframe
		replacementList[[i]] <- phylocom
	}
	
	#reduce the list to a single data frame and convert to matrix
	newCDM <- Reduce(rbind, replacementList)
	newCDM <- sample2matrix(newCDM)
	
	#it comes out in a weird order, so sort back to same plot order as input CDM
	newCDM <- newCDM[row.names(picante.cdm),]
	
	#add columns for species that were not recorded in any plots
	notFound <- setdiff(tree$tip.label, names(newCDM))
	
	if(length(notFound > 0))
	{
		toBind <- matrix(nrow=dim(newCDM)[[1]], ncol=length(notFound), 0)
		colnames(toBind) <- notFound
		newCDM <- cbind(newCDM, toBind)
	}
	
	#and sort into same species order as input cdm
	newCDM <- newCDM[,colnames(picante.cdm)]
	
	#for reasons that are not entirely clear to me, this is still returning as a DF
	#convert to matrix
	newCDM <- as.matrix(newCDM)
	
	newCDM
}

selectNear <- function(distances.between)
{
	#distances.between is a vector of distances between the focal cell and other cells
	#first exclude distances to the focal cell (and any other with distance = 0)
	distances.between <- distances.between[distances.between != 0]
	
	#now sample a cell name with a probability proportional to the inverse of the distance
	#from the focal cell
	newCell <- sample(x=names(distances.between), size=1, prob=1/distances.between)
	
	newCell
}

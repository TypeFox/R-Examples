#' Confirm that a CDM will run
#'
#' Check that a CDM will work with the dispersalNull model.
#'
#' @param picante.cdm A picante-style community data matrix with sites as rows, and
#' species as columns
#' 
#' @details It is possible that a CDM is unsuitable for the dispersalNull model. 
#' Specifically, if any single grid cell contains more species than are contained in the
#' sum of the remaining grid cells, the model will get stuck in an indefinite loop.
#'
#' @return either "pass" or "fail"
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
#' sim.abundances <- round(rlnorm(5000, meanlog=2, sdlog=1)) + 1
#'
#' cdm <- simulateComm(tree, richness.vector=10:25, abundances=sim.abundances)
#'
#' #this CDM should pass the check
#' checkCDM(cdm)
#'
#' #this CDM should not pass the check
#' test <- matrix(nrow=3, ncol=10)
#' test[1,] <- 1:10
#' test[2,] <- c(1,1,0,0,0,0,0,0,0,0)
#' test[3,] <- c(1,1,0,0,0,0,0,0,0,0)
#' test <- as.data.frame(test)
#' names(test)<-paste("s",1:10, sep="")
#' row.names(test) <- c("cell1","cell2","cell3")
#' checkCDM(test)

checkCDM <- function(picante.cdm)
{
	#set up a blank vector to save results into
	results <- c()

	for(i in 1:dim(picante.cdm)[1])
	{
		#for each row in the cdm, create a temporary CDM that is missing that row
		tempCDM <- picante.cdm[-i,]
		#calculate the observed species richness in the row in question
		obsRichness <- apply(picante.cdm, 1, lengthNonZeros)[i]
		#then find the total species richness in this temporary CDM. first sum all
		#occurrences, then find the length of non zeros
		tempRichness <- apply(tempCDM, 2, sum)
		otherRichness <- lengthNonZeros(tempRichness)
		#if the observed richness in the row in question is greater than otherRichness
		#set the relevant element of results to TRUE
		if(obsRichness > otherRichness)
		{
			results[i] <- TRUE
		}
		else
		{
			results[i] <- FALSE
		}
	}
	
	#test if any elements are TRUE in results and return either pass or fail
	if(any(results))
	{
		output <- "fail"
	}
	else
	{
		output <- "pass"
	}
	output
}

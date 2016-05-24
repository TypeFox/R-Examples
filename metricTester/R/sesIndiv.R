#' Summary statistics of SES results
#'
#' Summarizes individual iteration performance of each sim/null/metric combination across
#' all iterations
#'
#' @param raw.results A list of lists of lists of dataframes; the results of a call to
#' readIn.
#' @param concat.by Whether randomizations were concatenated by richness, plot or both.
#'
#' @details This function takes a raw list of results from multiple iterations from
#' multiLinker, and runs sesSingle across each one. It then summarizes the results of
#' those single runs as the number of sim/null/metrics that deviated beyond expectations
#' and the number that were within expectations. A single run from a given unique metric +
#' null approach is considered as throwing a type I error only if p is less
#' than or equal to 0.05 for the random spatial simulation. It would be possible to also
#' assess whether such unique combinations throw the opposite signal than expected for
#' habitat filtering and competitive exclusion. A unique combination iteration is
#' considered to throw a type II error if the p value from  either the filtering or the
#' exclusion simulation is greater than 0.05. 
#'
#' @return A data frame summarizing the total number of runs per spatial simulation, null,
#' metric, concat.by combination, as well as the type I and II error rates of each such
#' unique combination. 
#'
#' @export
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' #not run
#' #results <- readIn()
#' #summ <- sesIndiv(results, "both")

sesIndiv <- function(raw.results, concat.by)
{
	#lapply sesSingle across the list of results (one element per iteration)
	temp <- lapply(raw.results, sesSingle, concat.by)
	
	#reduce the list of results to one huge data frame
	temp <- Reduce(rbind, temp)
	
	#find the unique sim/null/metric combinations and turn into a temporary data frame
	#to loop over
	output <- unique(temp[,1:4])
	
	#set up an empty vector of Type I errors. as currently implemented, this will only
	#be triggered by runs where the random was considered significant. there are no TypeI
	#errors (currently) for filtering or competition. just TypeII. also set up an empty
	#vector to keep track of how many iterations there were for the given sim/null/metric
	
	totalRuns <- c()
	typeI <- c()
	typeII <- c()
	
	for(i in 1:dim(output)[1])
	{
		#this is a complicated subsetting operation. subset the reduced DF of results
		#(temp) down to instances where the sim, null, concat.by & metric match one
		#of the unique rows (i) from output
		temp2 <- temp[temp$simulation %in% output$simulation[i] & 
			temp$null.model %in% output$null[i] &
			temp$concat.by %in% output$concat.by[i] &
			temp$metric %in% output$metric[i], ]
		totalRuns[i] <- dim(temp2)[1]
		if(output$simulation[i] == "random")
		{
			typeII[i] <- NA
			typeI[i] <- dim(temp2[temp2$p.value <= 0.05,])[1]
		}
		else
		{
			typeI[i] <- NA
			typeII[i] <- dim(temp2[temp2$p.value > 0.05,])[1]
		}
	}
	
	#bind the three vectors into the output and send out
	output <- data.frame(output, total.runs=totalRuns, typeI, typeII)
	
	output
}

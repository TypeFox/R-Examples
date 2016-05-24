#' Read in the results of multiple metric/null/simulation tests
#'
#' This function is used to read all .RDS files from the working directory into a
#' single list. This list is then summarized with additional functions.
#'
#' @param working.directory Optional character string specifying the working directory.
#' If missing, the current working directory will be used. 
#' 
#' @details This function reads all .RDS files from the working directory (or another
#' specified directory) into a single list.
#'
#' @return A list of simulation results. Each element in the list consists of a single
#' output from multiLinker(). 
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

readIn <- function(working.directory)
{
	if(missing(working.directory))
	{
		working.directory <- getwd()
	}

	#return a character vector of all the files in the working directory
	allFiles <- list.files(path=working.directory)
	
	#subset to just those files that .RDS
	files <- allFiles[grep(".RDS", allFiles)]
	
	#paste the working directory into the name
	files <- paste(working.directory, files, sep="/")
	
	#save each rds file into a different element of the output list
	results <- lapply(files, readRDS)
	
	names(results) <- paste("iteration", 1:length(results), sep="")
	
	return(results)
}

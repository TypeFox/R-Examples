#' Overall per simulation-null-metric plot test
#'
#' This function provides one of many ways of summarizing and considering simulation
#' results.
#'
#' @param simulation.list A summarized results list such as one output from
#' reduceResults(). See examples.
#' @param concat.by Whether randomizations were concatenated by richness, plot or both.
#'
#' @details This function provides one way of summarizing and considering simulation
#' results. It takes as input a vector of 0s, 1s and 2s (corresponding to within, less,
#' and greater than the 95\% CIs, respectively) for all plots
#' from a given simulation-null-metric combination, and determines how many plots
#' overall deviated beyond expectations. A number of utility functions used for that goal
#' are also defined but not exported in this function. 
#'
#' @return A data frame summarizing the total number of plots tested and how many of
#' these deviated above (significantly overdispersed) or below (significantly clustered)
#' the 95\% CI for each simulation, null, metric combination. 
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
#' #summ <- reduceResults(results, "both")
#' #test <- plotOverall(summ$plot, "both")

plotOverall <- function(simulation.list, concat.by)
{
	if(concat.by=="richness" | concat.by=="plot")
	{
		#lapply all ones, twos and length LApply functions over the entire simulation list
		tempOnes <- lapply(simulation.list, onesLApply)
		tempTwos <- lapply(simulation.list, twosLApply)
		tempLength <- lapply(simulation.list, lengthLApply)
	
		#unlist and bind together into a data frame
		ones <- unlist(tempOnes)
		twos <- unlist(tempTwos)
		lengths <- unlist(tempLength)
	
		temp <- data.frame(ones, twos, lengths)
	
		#create better categories (i.e. separate the single string names into simulation,
		#null and metric columns). this is a long command but it runs really quickly
		betterNames <- suppressWarnings(data.frame(Reduce(rbind, 
			strsplit(names(ones), split="[.]"))))
		betterNames$concat.by <- concat.by
		
		output <- cbind(betterNames, temp)

		names(output) <- c("simulation", "null.model", "metric", "concat.by",
			"clustered", "overdispersed", "total.plots")
	
		row.names(output) <- NULL
	}

	else if(concat.by=="both")
	{
		#elaborate nested lapply statement to run lengths LApply functions over entire
		#summarized results
		ones <- lapply(seq_along(simulation.list), function(x)
			unlist(lapply(simulation.list[[x]], onesLApply)))
		twos <- lapply(seq_along(simulation.list), function(x)
			unlist(lapply(simulation.list[[x]], twosLApply)))
		lengths <- lapply(seq_along(simulation.list), function(x)
			unlist(lapply(simulation.list[[x]], lengthLApply)))
		#this did not give names to elements of the resulting list, so give them names
		#so that when you unlist those names will get carried along
		names(ones) <- names(simulation.list)
		names(twos) <- names(simulation.list)
		names(lengths) <- names(simulation.list)
		ones <- unlist(ones)
		twos <- unlist(twos)
		lengths <- unlist(lengths)
		#now come up with better names, like you did when concat was richness or plot
		#however, the by.richness or q will mess things up, so will need to delete those
		betterNames <- suppressWarnings(data.frame(Reduce(rbind, 
			strsplit(names(ones), split="[.]"))))
		betterNames[,3] <- NULL
		
		output <- data.frame(betterNames, ones, twos, lengths)

		names(output) <- c("simulation", "null.model", "concat.by", "metric",
			"clustered", "overdispersed", "total.plots")
		row.names(output) <- NULL
	}

	else
	{
		stop("concat.by must equal either both, richness, or plot")
	}
	
	output
}

lengthApply <- function(dataframe)
{
	#this utility function will tell you how many of each simulation got run
	#exclude "richness" and "plot" columns
	exclude <- c("richness", "plot")
	temp <- dataframe[ ,!(names(dataframe) %in% exclude)]

	#apply length
	output <- apply(temp, 2, length)
	
	output
}

lengthLApply <- function(null.list)
{
	#lapply tWrapApply over null.list
	temp <- lapply(null.list, lengthApply)

	#unlist the output into a single vector. let it assign names, it does a decent job
	output <- unlist(temp)

	#just return the simple vector with ugly names, fix the names in a higher function
	#later
	output
}

lengthOnes <- function(input.vector)
{
	#ones are significantly clustered
	ones <- input.vector[input.vector == 1]
	return(length(ones))
}

lengthTwos <- function(input.vector)
{
	#twos are significantly overdisperesed
	twos <- input.vector[input.vector == 2]
	return(length(twos))
}

onesApply <- function(dataframe)
{
	#exclude "richness" and "plot" columns
	exclude <- c("richness", "plot")
	temp <- dataframe[ ,!(names(dataframe) %in% exclude)]

	#apply tWrap over a data frame of metric SES scores for a given null and spatial sim
	output <- apply(temp, 2, lengthOnes)
	
	output
}

onesLApply <- function(null.list)
{
	#lapply tWrapApply over null.list
	temp <- lapply(null.list, onesApply)

	#unlist the output into a single vector. let it assign names, it does a decent job
	output <- unlist(temp)

	#just return the simple vector with ugly names, fix the names in a higher function
	#later
	output
}

twosApply <- function(dataframe)
{
	#exclude "richness" and "plot" columns
	exclude <- c("richness", "plot")
	temp <- dataframe[ ,!(names(dataframe) %in% exclude)]

	#apply lengthTwos over a df of metric SES scores for a given null and spatial sim
	output <- apply(temp, 2, lengthTwos)
	
	output
}

twosLApply <- function(null.list)
{
	#lapply tWrapApply over null.list
	temp <- lapply(null.list, twosApply)

	#unlist the output into a single vector. let it assign names, it does a decent job
	output <- unlist(temp)

	#just return the simple vector with ugly names, fix the names in a higher function
	#later
	output
}

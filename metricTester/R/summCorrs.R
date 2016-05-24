#' Summarize correlations among metrics over a result from a varyX function
#'
#' Takes the results of one of the varyX functions, and calculates the correlations among
#' metrics, returning either the raw or summarized correlations.
#'
#' @param vary.results Results from any of the varyX (e.g., varyAbundance) functions.
#' @param exclude Results columns to exclude from correlation. For instance, with alpha
#' metrics, one would want to, at the minimum, exclude the plot name column.
#' @param return.raw Default is FALSE. Whether to return the raw correlation coefficients
#' between the metrics for each element from vary.results, or whether to summarize the
#' correlations by their mean per parameter set from vary.results.
#' @param cor.method Default is "spearman", but takes any of the other options from the
#' base cor function.
#' 
#' @details Not a well tested function.
#'
#' @return If return.raw is set to FALSE, returns a list of matrices, one for each set
#' of parameters in vary.results, summarizing the mean correlation coefficients between
#' each pairwise metric correlation. If return.raw is set to TRUE, returns a list of
#' lists, one for each set of parameters. Each of these lists is the length of the number
#' of iterations in vary.results. Each element in the list is a matrix of pairwise metric
#' correlations. 
#'
#' @export
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' #below not run for timing issues on CRAN
#' #system.time(vSize <- varyX(alpha=TRUE, tree.size=c(40, 50), richness=20:30, delta=1,
#'	#abundances=round(rlnorm(5000, meanlog=2, sdlog=1)) + 1, iterations=2, cores="seq"))
#'
#' #test <- summCorrs(vSize, exclude=c("plot", "richness"))

summCorrs <- function(vary.results, exclude, return.raw=FALSE, cor.method="spearman")
{
	metrics <- names(vary.results[[1]][[1]])
	
	metrics <- metrics[!(metrics %in% exclude)]

	allResults <- list()
	#i refers to iterations
	for(i in 1:length(vary.results))
	{
		#this will get reset each iteration
		singleSet <- list()
		
		#j refers to a single one of the set of varied options (e.g. a single tree size)
		for(j in 1:length(vary.results[[i]]))
		{
			singleSet[[j]] <- suppressWarnings(cor(vary.results[[i]][[j]][,metrics],
				method=cor.method))
		}
		
		allResults[[i]] <- singleSet
	}
	
	#make a list that will end up being the length of the number of varied options
	#so, i here refers to a varied option
	betterOrganized <- list()
	
	for(i in 1:length(vary.results[[1]]))
	{
		#make a list here equal to the number of iterations
		iterationList <- list()
		
		#j refers to iterations
		for(j in 1:length(vary.results))
		{
			iterationList[[j]] <- allResults[[j]][[i]]
		}
	
		names(iterationList) <- names(vary.results)
		betterOrganized[[i]] <- iterationList
	}
	
	names(betterOrganized) <- names(vary.results[[1]])
	
	if(return.raw==FALSE)
	{
		results <- list()
		
		for(i in 1:length(betterOrganized))
		{
			singleCorr <- apply(simplify2array(betterOrganized[[i]]), 1:2, mean, na.rm=T)
			results[[i]] <- singleCorr
		}

		names(results) <- names(vary.results[[1]])
	}
	
	else if(return.raw==TRUE)
	{
		results <- betterOrganized
	}

	else
	{
		stop("return.raw must be either TRUE or FALSE")
	}
	
	results
}

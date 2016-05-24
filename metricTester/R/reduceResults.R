#' Reduce results from multiLinker into a manageable format
#'
#' Calling multiLinker creates .RDS files, one per iteration. This function will
#' combine these results into a more manageable format.
#'
#' @param results.list The results of a call to readIn()
#' @param concat.by Whether randomizations were concatenated by richness, plot or both
#'
#' @details Given a list of results readIn() from multiLinker, this function will reduce
#' those results into a manageable format like that expected for calls to plotOverall
#' and sesOverall.
#'
#' @return A list of data frames. 
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

reduceResults <- function(results.list, concat.by)
{
	if(!(concat.by %in% c("both","plot","richness")))
	{
		stop("concat.by must equal either both, richness, or plot")
	}

	#assume that all iterations have same dimensions. should ultimately write a check here
	#pull the iteration, sim, null and metric names out for use later
	iterations <- names(results.list)
	sims <- names(results.list[[1]])
	nulls <- names(results.list[[1]][[1]][[1]])

	#combine each element from the long list of results, where each is one iteration, 
	#into a list where each element from the first level of the list corresponds to a
	#spatial simulation. each second level are the combined results of a single iteration
	#from multiLinker. this second level is a matrix of lists, where the first
	#column relates to the ses, the second to the plot significance. there are as many
	#rows in the matrix as there are iterations from multiLinker
	firstLevel <- reduceRandomizations(results.list)
	
	#pull the arena and plotTest results out separately
	ses <- list()
	
	for(i in 1:length(firstLevel))
	{
		ses[[i]] <- firstLevel[[i]][,1]
	}
	
	plot <- list()
	
	for(i in 1:length(firstLevel))
	{
		plot[[i]] <- firstLevel[[i]][,2]
	}
	
	#give those separate results names for the simulations
	names(ses) <- sims
	names(plot) <- sims

	secondLevel <- list("ses"=ses, "plot"=plot)
	
	#if concat.by=both, just brute force the results into an acceptable order
	if(concat.by=="both")
	{
		#set up an empty list 2 long, one for ses one for plot
		results <- vector("list", 2)
		names(results) <- c("ses", "plot")
		#i elements refer to either ses or plot
		for(i in 1:length(secondLevel))
		{
			#set up an empty list as long as the number of spatial sims
			jTemp <- vector("list", length(sims))
			names(jTemp) <- sims
			#j elements refer to spatial simulations
			for(j in 1:length(sims))
			{
				#k elements refer to the number of null models
				kTemp <- vector("list", length(nulls))
				names(kTemp) <- nulls
				for(k in 1:length(nulls))
				{
					#set up two empty lists to pull each iteration of richness and plot
					#from a given null in as a data frame
					richnessTemp <- vector("list", length(nulls))
					plotTemp <- vector("list", length(nulls))
					
					#l elements refer to the number of iterations
					for(l in 1:length(iterations))
					{
						#pull all the richness and plot data frames for a given null
						#out and make each into a new element in a list
						richnessTemp[[k]][[l]] <- secondLevel[[i]][[j]][[l]][[k]]$richness
						plotTemp[[k]][[l]] <- secondLevel[[i]][[j]][[l]][[k]]$plot
					}
					richnessReduced <- Reduce(rbind, richnessTemp[[k]])
					plotReduced <- Reduce(rbind, plotTemp[[k]])
					kTemp[[k]] <- list("by.richness"=richnessReduced,
						"by.plot"=plotReduced)
				}
			jTemp[[j]] <- kTemp
			}
		results[[i]] <- jTemp
		}
	}
	
	else
	{
		results <- lapply(secondLevel, function(x) lapply(x, reduceRandomizations))
	}
	
	results
}

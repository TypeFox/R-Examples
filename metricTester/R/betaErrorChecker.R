#' Wrapper for summarizing error rates of beta metric randomizations
#'
#' Given the results of a single iteration of the betaLinker function, returns a list of data
#' frames summarizing the type I and II error rates of 
#' metrics both at the single plot and the entire arena level.
#'
#' @param single.iteration Results of a run of the betaLinker function.
#' 
#' @details This function wraps a number of smaller functions into a useful type I and II
#' error checker. It takes a reduced list of randomizations such as those reduced from
#' metricsNnulls with reduceRandomizations, summarizes the mean,
#' SD, and CI of each metric plus null model either at the richness or plot level,
#' then compares the observed metric scores to those summarized metrics. It return a list
#' with two elements. The first is a list of data frames, where each corresponds to the 
#' standardized effect scores of the observed metrics for a given null model. The second
#' is a list of data frames, where each corresponds to whether a given plot deviates
#' beyond CI. For the latter, 0 corresponds to within CI, 1 corresponds to less than the
#' CI, and 2 corresponds to greater than the CI.
#'
#' @return A list of lists of matrices. The first element of the lists refers to the
#' results from a spatial simulation. Within each of these elements is a list of matrices,
#' where each matrix tabulates the error rate of all tested beta diversity metrics with
#' a given null model.
#'
#' @export
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' #run the betaLinker function
#' #below not run for timing issues on CRAN
#' #system.time(ex <- betaLinker(no.taxa=50, arena.length=300, mean.log.individuals=2, 
#' #length.parameter=5000, sd.parameter=50, max.distance=30, proportion.killed=0.2, 
#'	#competition.iterations=3, no.plots=15, plot.length=30,
#'	#randomizations=3, cores="seq",
#'	#nulls=list("richness"=metricTester:::my_richnessNull,
#'	#"frequency"=metricTester:::my_frequency)))
#'
#' #test <- betaErrorChecker(ex)

betaErrorChecker <- function(single.iteration)
{
	#this is a cheap hack to figure out which metrics (and below which nulls) got used.
	#we need to ensure that all metrics used have an expectation. for instance, it
	#is possible to define your own beta metric and use that, but it would not have an
	#expectation here. Metrics like Ist signify clustering when the observed value is
	#greater than the randomized values, while metrics like MPD signify overdispersion
	#when the observed value is greater than the randomized values. to be able to use the
	#same greater or less than code below, flip the expectations for Ist and any other
	#metrics where greater than the randomized values equals filtering.
	
	usedMetrics <- names(single.iteration$observed[[1]])
	#pull out richness, total abundance, and any other functions (like skewness)
	#that are not phylo measures
	usedMetrics <- usedMetrics[usedMetrics != "richness" & usedMetrics!="total_abundance"]
	
	reversed <- c("Ist", "Pst", "Bst", "PIst", "Ist_unpruned", "Pst_unpruned",
		"Bst_unpruned", "PIst_unpruned")
	
	standard <- usedMetrics[!(usedMetrics %in% reversed)]
	
	usedNulls <- names(single.iteration$randomized[[1]])

	#set up a blank matrix to save results in. will re-use this matrix for each spatial
	#simulation
	resultMatrix <- matrix(nrow=1, ncol=length(usedMetrics))

	#set up a blank list to save results of blank matrix into as it gets filled. note that
	#because we have the option to set elements in result matrix to 0, we do not need to
	#reset it after each iteration. to facilitate moving down the list with each element
	#of k, just start a placeholder here
	results <- list()
	placeholder <- 0

	#i refers to simulations. set up a sublist for each spatial simulation, where each
	#element of that list refers to a given null model
	for(i in 1:length(single.iteration$observed))
	{
		temp <- list()
		#j refers to null models.
		for(j in 1:length(usedNulls))
		{
			#k refers to metrics
			for(k in 1:length(usedMetrics))
			{
				if(names(single.iteration$observed)[i]=="random")
				{
					if(single.iteration$observed[[i]][usedMetrics[k]] < 
						quantile(single.iteration$randomized[[i]][[usedNulls[j]]][,usedMetrics[k]], 0.025))
					{
						resultMatrix[1, k] <- 1
					}
					else if(single.iteration$observed[[i]][usedMetrics[k]] > 
						quantile(single.iteration$randomized[[i]][[usedNulls[j]]][,usedMetrics[k]], 0.975))
					{
						resultMatrix[1, k] <- 1
					}
					else
					{
						resultMatrix[1, k] <- 0
					}
					rownames(resultMatrix) <- "typeI"
				}
				else if(names(single.iteration$observed)[i]=="filtering" & usedMetrics[k] %in% reversed)
				{
					if(single.iteration$observed[[i]][usedMetrics[k]] <
						quantile(single.iteration$randomized[[i]][[usedNulls[j]]][,usedMetrics[k]], 0.95))
					{
						resultMatrix[1, k] <- 1
					}
					else
					{
						resultMatrix[1, k] <- 0
					}
					rownames(resultMatrix) <- "typeII"
				}
				else if(names(single.iteration$observed)[i]=="filtering" & usedMetrics[k] %in% standard)
				{
					if(single.iteration$observed[[i]][usedMetrics[k]] >
						quantile(single.iteration$randomized[[i]][[usedNulls[j]]][,usedMetrics[k]], 0.05))
					{
						resultMatrix[1, k] <- 1
					}
					else
					{
						resultMatrix[1, k] <- 0
					}
					rownames(resultMatrix) <- "typeII"
				}
				else if(names(single.iteration$observed)[i]=="competition" & usedMetrics[k] %in% reversed)
				{
					if(single.iteration$observed[[i]][usedMetrics[k]] >
						quantile(single.iteration$randomized[[i]][[usedNulls[j]]][,usedMetrics[k]], 0.05))
					{
						resultMatrix[1, k] <- 1
					}
					else
					{
						resultMatrix[1, k] <- 0
					}
					rownames(resultMatrix) <- "typeII"
				}
				else if(names(single.iteration$observed)[i]=="competition" & usedMetrics[k] %in% standard)
				{
					if(single.iteration$observed[[i]][usedMetrics[k]] <
						quantile(single.iteration$randomized[[i]][[usedNulls[j]]][,usedMetrics[k]], 0.95))
					{
						resultMatrix[1, k] <- 1
					}
					else
					{
						resultMatrix[1, k] <- 0
					}
					rownames(resultMatrix) <- "typeII"
				}
			}
		placeholder <- placeholder + 1
		colnames(resultMatrix) <- usedMetrics
		temp[[j]] <- resultMatrix
		}
		names(temp) <- usedNulls
		results[[i]] <- temp
	}
	
	names(results) <- names(single.iteration$observed)
	results
}

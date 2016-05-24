#' Overall per simulation-null-metric SES test
#'
#' This function provides one of many ways of summarizing and considering simulation
#' results.
#'
#' @param simulation.list A summarized results list such as one output from
#' reduceResults(). See examples.
#' @param test Either "ttest" or "wilcotest", depending on whether the user wants to run
#' a two-sided t-test or a Wilcoxon signed rank test.
#' @param concat.by Whether randomizations were concatenated by richness, plot or both.
#'
#' @details This function provides one way of summarizing and considering simulation
#' results. It takes as input a vector of all standardized effect sizes for all plots
#' from a given simulation-null-metric combination, and calculates the mean of the vector
#' and whether it differs significantly from a mean of zero. It does this either with a
#' simple two-sided t-test, or with a Wilcoxon signed rank test. If the latter, and if
#' there are three different spatial simulations with names random, filtering and
#' competition, the test is two-sided, less and greater, respectively. 
#'
#' @return A data frame summarizing the mean, overall standardized effect sizes and the
#' significance of those devations from expectations for each simulation, null, metric
#' combination. This test works across all iterations, and looks for overall shifts in
#' SES from expectations (see details for for expectations).
#'
#' @export
#'
#' @importFrom dplyr select
#' @importFrom stats cor dist sd t.test weighted.mean wilcox.test
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' #not run
#' #results <- readIn()
#' #summ <- reduceResults(results)
#' #examp <- sesOverall(summ$ses, "both")

sesOverall <- function(simulation.list, test, concat.by)
{
	#dumb hack to pass R CMD check
	simulation <- "hack"
	null.model <- "hack"
	metric <- "hack"
	estimate <- "hack"
	p.value <- "hack"

	if(!(concat.by %in% c("both","plot","richness")))
	{
		stop("concat.by must equal either both, richness, or plot")
	}

	#add a line to throw an error if test is not properly specified
	if(test != "ttest" & test != "wilcotest")
	{
		stop("test must be set to one of 'ttest' or 'wilcotest'")
	}

	#if test is set to ttest, just do a simple t.test to see if mean differs from 0.
	if(test=="ttest" & (concat.by=="richness" | concat.by=="plot"))
	{
		#lapply tWrapLApply over simulation.list
		tempAll <- lapply(simulation.list, tWrapLApply)
	}
	
	else if(test=="ttest" & concat.by=="both")
	{
		stop("t.test with concat by both not currently operational")
	}
	
	#more probably people will use the wilcoxon signed rank test. go into that here
	else if(test=="wilcotest" & (concat.by=="richness" | concat.by=="plot"))
	{
		#if the names of the spatial simulations are "random" "competition" and 
		#"filtering", create a character vector of expected alternative hypotheses to feed
		#into anonymous function below
		if(setequal(names(simulation.list), c("random","filtering","competition")))
		{
			toFeed <- names(simulation.list)
			toFeed[toFeed=="random"] <- "two.sided"
			toFeed[toFeed=="filtering"] <- "less"
			toFeed[toFeed=="competition"] <- "greater"
		}
		
		else if(!setequal(names(simulation.list), c("random","filtering","competition")))
		{
			print("You included new spatial simulations. Modify expectations manually")
			toFeed <- rep("two.sided", 3)
		}
		
		tempAll <- lapply(seq_along(toFeed), function(x)
			wilcoWrapLApply(simulation.list[[x]], alternative=toFeed[x]))
	}
	
	else if(test=="wilcotest" & concat.by=="both")
	{
		if(setequal(names(simulation.list), c("random","filtering","competition")))
		{
			toFeed <- names(simulation.list)
			toFeed[toFeed=="random"] <- "two.sided"
			toFeed[toFeed=="filtering"] <- "less"
			toFeed[toFeed=="competition"] <- "greater"
		}
		
		else if(!setequal(names(simulation.list), c("random","filtering","competition")))
		{
			print("You included new spatial simulations. Modify expectations manually")
			toFeed <- rep("two.sided", 3)
		}

		tempAll <- list()
		
		#i refers to spatial simulation (you will only feed in $ses components)
		for(i in 1:length(toFeed))
		{
			#j refers to null models
			tempNull <- list()
			for(j in 1:length(simulation.list[[1]]))
			{
				#k refers either to concatenating by richness or by plot
				for(k in 1:length(simulation.list[[1]][[1]]))
				{
					temp <- wilcoWrapLApply(simulation.list[[i]][[j]],
						alternative=toFeed[i])
					#this function was designed for use over a list of null models, so
					#change name here to reflect what it is actually being used for
					names(temp)[4] <- "concat.by"
				}
				#set the jth element of tempNull equal to the dataframe temp
				tempNull[[j]] <- temp
				#repeat the name of the null the correct length of the null model and bind
				tempNull[[j]]$null.model <- rep(names(simulation.list[[i]][j]),
					dim(tempNull[[j]])[1])
			}
			#set the ith element of tempAll equal to the list tempNull
			tempAll[[i]] <- tempNull
			#reduce the list down
			tempAll[[i]] <- Reduce(rbind, tempAll[[i]])
			#repeat the sim name the correct length and bind
			tempAll[[i]]$simulation <- rep(names(simulation.list[i]),
				dim(tempAll[[i]])[1])
		}
	}
	
	#reduce the output list into a single data frame
	output <- Reduce(rbind, tempAll)
	
	if(concat.by=="richness" | concat.by=="plot")
	{
		#create a vector of expanded simulation names. note that this code is sensitive to
		#changes. for instance, if one simulation tests certain nulls that another does
		#not, this will not end up being correct. this generates a data frame, but we only
		#save the second column
		simNames <- expand.grid(tempAll[[1]]$null.model, names(simulation.list))[,2]
	
		output$simulation <- simNames
	
		output <- cbind(output, concat.by=rep(concat.by, dim(output)[1]))
	
		output <- select(output, simulation, null.model, metric, concat.by,
			estimate, p.value)		
	}
	
	else if(concat.by=="both")
	{
		output <- select(output, simulation, null.model, metric, concat.by,
			estimate, p.value)
	}
	
	output
}

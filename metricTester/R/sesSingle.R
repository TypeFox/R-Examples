#' Summary of results from a single iteration
#'
#' Use Wilcoxon signed rank test to determine whether plots from a SINGLE iteration 
#' differ from expectations
#'
#' @param single.iteration The results of a single iteration from multiLinker.
#' @param concat.by Whether randomizations were concatenated by richness, plot or both.
#'
#' @details This function uses a Wilcoxon signed rank test to determine whether the
#' plots from a spatial simulation/null/metric from a SINGLE iteration differ from
#' expectations. Assuming there are three spatial simulations named random, filtering, and
#' competition, this function will use two.sided, lesser and greater Wilcoxon tests,
#' respectively.
#'
#' @return A data frame summarizing the mean of standardized effect sizes and the
#' significance of those devations from expectations for a given iteration (e.g. the 
#' plots from a given arena). It does consider all spatial simulations, nulls and 
#' metrics from that iteration. This test works across a single iteration, and will
#' generally not be used by itself; it is called by sesIndiv.
#'
#' @export
#'
#' @importFrom dplyr select
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' #not run
#' #results <- readIn()
#' #summ <- sesSingle(results[[1]], "richness")

sesSingle <- function(single.iteration, concat.by)
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

	#identify any null model/spatial sim combos that failed. 
	problems <- failed(single.iteration, concat.by)
	
	#if there were any problems use a for loop to go through each failed run and set that 
	#element to NULL in the single.iteration results
	if(dim(problems)[1] > 0)
	{
		for(i in 1:dim(problems)[1])
		{
			sim <- problems[i,"simulation"]
			null <- problems[i,"null"]
			if(concat.by=="richness" | concat.by=="plot")
			{
				single.iteration[[sim]]$ses[[null]] <- NULL
			}
			else
			{
				single.iteration[[sim]]$ses[[null]]$richness <- NULL
			}
		}
	}

	#if the names of the spatial simulations are "random" "competition" and 
	#"filtering", create a character vector of expected alternative hypotheses to feed
	#into anonymous function below
	if(setequal(names(single.iteration), c("random","filtering","competition")))
	{
		toFeed <- names(single.iteration)
		toFeed[toFeed=="random"] <- "two.sided"
		toFeed[toFeed=="filtering"] <- "less"
		toFeed[toFeed=="competition"] <- "greater"
	}
		
	else if(!setequal(names(single.iteration), c("random","filtering","competition")))
	{
		print("You included new spatial simulations. Modify expectations manually")
		toFeed <- rep("two.sided", 3)
	}
	
	if(concat.by=="richness" | concat.by=="plot")
	{
		temp <- lapply(seq_along(toFeed), function(x)
			wilcoWrapLApply(single.iteration[[x]]$ses, alternative=toFeed[x]))

		#reduce the output list into a single data frame
		output <- Reduce(rbind, temp)
	
		#create a data frame of expanded simulation names
		simNames <- expand.grid(temp[[1]]$null.model, names(single.iteration))
	
		#give it names
		names(simNames) <- c("null", "simulation")
	
		#if there were any problems, take these names out of simNames
		if(dim(problems)[1] > 0)
		{
			simNames <- simNames[!(simNames$null %in% problems$null & 
				simNames$simulation %in% problems$simulation),]
		}
	
		#use just the second col from simNames
		output$simulation <- simNames[,2]
		
		#tack on concat names
		output$concat.by <- concat.by
	}
	
	else
	{
		tempAll <- list()
		
		#i refers to spatial simulations (you will only feed in $ses components)
		for(i in 1:length(toFeed))
		{
			#j refers to null models
			tempNull <- list()
			for(j in 1:length(single.iteration[[1]]$ses))
			{
				#k refers either to concatenating by richness or by plot
				for(k in 1:length(single.iteration[[1]]$ses[[1]]))
				{
					temp <- wilcoWrapLApply(single.iteration[[i]]$ses[[j]],
						alternative=toFeed[i])
					#this function was designed for use over a list of null models, so
					#change name here to reflect what it is actually being used for
					names(temp)[4] <- "concat.by"
				}
				#set the jth element of tempNull equal to the dataframe temp
				tempNull[[j]] <- temp
				#repeat the name of the null the correct length of the null model and bind
				tempNull[[j]]$null.model <- rep(names(single.iteration[[1]]$ses)[j],
					dim(tempNull[[j]])[1])
			}
			#set the ith element of tempAll equal to the list tempNull
			tempAll[[i]] <- tempNull
			#reduce the list down
			tempAll[[i]] <- Reduce(rbind, tempAll[[i]])
			#repeat the sim name the correct length and bind
			tempAll[[i]]$simulation <- rep(names(single.iteration)[i],
				dim(tempAll[[i]])[1])
		}
		#reduce the output list into a single data frame
		output <- Reduce(rbind, tempAll)	
	}
	
	#select the columns so they come out in a nice order
	output <- select(output, simulation, null.model, concat.by, metric,
		estimate, p.value)
	
	output
}

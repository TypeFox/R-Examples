#' Return summary statistics from a data frame of randomized metric values
#'
#' Summarizes observed metric scores. Returns the mean, standard deviation and 
#' 95\% confidence intervals of each plot or observed richness. 
#'
#' @param null.output Data frame of randomized metric values such as an element from a
#' call to reduceRandomizations()
#' @param concat.by Whether to concatenate the randomizations by richness, plot or both
#' 
#' @details Given a data frame of metric values, summarizes either by plot or richness.
#' Outputs the mean, standard deviation and 95\% confidence intervals of each plot or
#' observed richness. If provided with concat.by="both", outputs a list of two data
#' frames, one for by richness and one for by plot. Otherwise, outputs a data frame.
#'
#' @return Either a list of or a data frame of summarized metric scores, see details.
#'
#' @export
#'
#' @importFrom dplyr group_by summarize
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
#' #below not run for example timing issues on CRAN
#'
#' #rawResults <- metricsNnulls(tree, cdm)
#'
#' #results <- reduceRandomizations(rawResults)
#'
#' #test <- summaries(results$frequency, concat.by="richness")

summaries <- function(null.output, concat.by="richness")
{
	#dumb hack to pass R CMD check
	richness <- "hack"
	metric <- "hack"

	#set up the results table to be the appropriate size if concat.by = richness
	#you want the table to have one row for every richness value and four columns for each
	#metric plus an additional column for the richness
	if(concat.by=="richness")
	{
		results <- matrix(ncol=4 * (dim(null.output)[2]-2) + 1, 
			nrow=length(unique(null.output$richness)), 0)
	}
	#set up the results table to be the appropriate size if concat.by = plot
	else if(concat.by=="plot")
	{
		results <- matrix(ncol=4 * (dim(null.output)[2]-2) + 1, 
			nrow=length(unique(null.output$plot)), 0)
	}
	else if(concat.by=="both")
	{
		results <- list()
		results$richness <- matrix(ncol=4 * (dim(null.output)[2]-2) + 1, 
			nrow=length(unique(null.output$richness)), 0)
		results$plot <- matrix(ncol=4 * (dim(null.output)[2]-2) + 1, 
			nrow=length(unique(null.output$plot)), 0)
	}
	#throw an error if none of above options
	else
	{
		stop("Argument concat.by must be 'richness', 'plot', or 'both'")
	}
	
	if(concat.by=="richness")
	{
		#start the for loop at i=3 to skip the richness and plot columns
		for(i in 3:dim(null.output)[2])
		{
			#create a temporary data frame because difficult to use dplyr over big tables
			temp <- data.frame(richness=null.output$richness, metric=null.output[,i])
			grouped <- group_by(temp, richness)
			#note that we want to start the last column at 5, to leave the first blank for
			#either plot or richness names
			lastCol <- (4*(i-2))+1
			results[,lastCol-3] <- as.data.frame(summarize(grouped, 
				mean(metric, na.rm = TRUE)))[,2]
			results[,lastCol-2] <- as.data.frame(summarize(grouped, 
				sd(metric, na.rm = TRUE)))[,2]
			results[,lastCol-1] <- as.data.frame(summarize(grouped, 
				quantile(metric, 0.025, na.rm=TRUE)))[,2]
			results[,lastCol] <- as.data.frame(summarize(grouped, 
				quantile(metric, 0.975, na.rm=TRUE)))[,2]
		}
		#dplyr automatically creates title for summarized columns, so just pull last used
		rows <- summarize(grouped, mean(metric, na.rm = TRUE))[,1]
	}
	else if(concat.by=="plot")
	{
		#start the for loop at i=3 to skip the richness and plot columns
		for(i in 3:dim(null.output)[2])
		{
			#create a temporary data frame because difficult to use dplyr over big tables
			temp <- data.frame(plot=null.output$plot, metric=null.output[,i])
			grouped <- group_by(temp, plot)
			#note that we want to start the last column at 5, to leave the first blank for
			#either plot or richness names
			lastCol <- (4*(i-2))+1
			results[,lastCol-3] <- as.data.frame(summarize(grouped, 
				mean(metric, na.rm = TRUE)))[,2]
			results[,lastCol-2] <- as.data.frame(summarize(grouped, 
				sd(metric, na.rm = TRUE)))[,2]
			results[,lastCol-1] <- as.data.frame(summarize(grouped, 
				quantile(metric, 0.025, na.rm=TRUE)))[,2]
			results[,lastCol] <- as.data.frame(summarize(grouped, 
				quantile(metric, 0.975, na.rm=TRUE)))[,2]
		}
		#dplyr automatically creates title for summarized columns, so just pull last used
		rows <- summarize(grouped, mean(metric, na.rm = TRUE))[,1]
	}
	else if(concat.by=="both")
	{
		#run the richness loop first
		for(i in 3:dim(null.output)[2])
		{
			#create a temporary data frame because difficult to use dplyr over big tables
			temp <- data.frame(richness=null.output$richness, metric=null.output[,i])
			grouped <- group_by(temp, richness)
			#note that we want to start the last column at 5, to leave the first blank for
			#either plot or richness names
			lastCol <- (4*(i-2))+1
			results$richness[,lastCol-3] <- as.data.frame(summarize(grouped, 
				mean(metric, na.rm = TRUE)))[,2]
			results$richness[,lastCol-2] <- as.data.frame(summarize(grouped, 
				sd(metric, na.rm = TRUE)))[,2]
			results$richness[,lastCol-1] <- as.data.frame(summarize(grouped, 
				quantile(metric, 0.025, na.rm=TRUE)))[,2]
			results$richness[,lastCol] <- as.data.frame(summarize(grouped, 
				quantile(metric, 0.975, na.rm=TRUE)))[,2]
		}
		#dplyr automatically creates title for summarized columns, so just pull last used
		rowsRichness <- summarize(grouped, mean(metric, na.rm = TRUE))[,1]
		
		#run the plot loop now
		#start the for loop at i=3 to skip the richness and plot columns
		for(i in 3:dim(null.output)[2])
		{
			#create a temporary data frame because difficult to use dplyr over big tables
			temp <- data.frame(plot=null.output$plot, metric=null.output[,i])
			grouped <- group_by(temp, plot)
			#note that we want to start the last column at 5, to leave the first blank for
			#either plot or richness names
			lastCol <- (4*(i-2))+1
			results$plot[,lastCol-3] <- as.data.frame(summarize(grouped, 
				mean(metric, na.rm = TRUE)))[,2]
			results$plot[,lastCol-2] <- as.data.frame(summarize(grouped, 
				sd(metric, na.rm = TRUE)))[,2]
			results$plot[,lastCol-1] <- as.data.frame(summarize(grouped, 
				quantile(metric, 0.025, na.rm=TRUE)))[,2]
			results$plot[,lastCol] <- as.data.frame(summarize(grouped, 
				quantile(metric, 0.975, na.rm=TRUE)))[,2]
		}
		#dplyr automatically creates title for summarized columns, so just pull last used
		rowsPlot <- summarize(grouped, mean(metric, na.rm = TRUE))[,1]
	}

	#add in the first column info
	if(concat.by == "richness" | concat.by== "plot")
	{
		#set the results to dataframe so it can have a column with plot names if needed
		results <- as.data.frame(results)
		results[,1] <- rows
	}
	else if(concat.by == "both")
	{
		#set the results to dataframe so it can have a column with plot names if needed
		results$richness <- as.data.frame(results$richness)
		results$plot <- as.data.frame(results$plot)
		results$richness[,1] <- rowsRichness
		results$plot[,1] <- rowsPlot
	}
	#create column names
	metricNames <- names(null.output)[names(null.output)!="richness" & 
		names(null.output)!="plot"]
	summaryNames <- c("average", "sd", "lower", "upper")
	#this last set of names is the metric name plus each of the different summaries
	comboNames <- paste(rep(metricNames, each = length(summaryNames)), 
		rep(summaryNames, length(metricNames)), sep = ".")
	#bind this to either "richness" or "plot" for name of first column
	if(concat.by=="richness")
	{
		comboNames <- c("richness", comboNames)
		names(results) <- comboNames
	}
	else if(concat.by=="plot")
	{
		comboNames <- c("plot", comboNames)
		names(results) <- comboNames
	}
	else if(concat.by=="both")
	{
		comboNamesRichness <- c("richness", comboNames)
		comboNamesPlot <- c("plot", comboNames)
		names(results$richness) <- comboNamesRichness
		names(results$plot) <- comboNamesPlot
	}
	return(results)
}

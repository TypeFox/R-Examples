#most non-exported utility functions are stored here

##########################################################################################
################################## ALPHA METRICS #########################################
##########################################################################################

my_richness <- function(metrics.input)
	apply(metrics.input$picante.cdm, 1, lengthNonZeros)

naw_mpd <- function(metrics.input)
	modifiedMPD(metrics.input$picante.cdm, metrics.input$dists, abundance.weighted=FALSE)

inter_mpd <- function(metrics.input)
	modifiedMPD(metrics.input$picante.cdm, metrics.input$dists, abundance.weighted="interspecific")

intra_mpd <- function(metrics.input)
	modifiedMPD(metrics.input$picante.cdm, metrics.input$dists, abundance.weighted="intraspecific")

complete_mpd <- function(metrics.input)
	modifiedMPD(metrics.input$picante.cdm, metrics.input$dists, abundance.weighted="complete")

naw_mntd <- function(metrics.input)
	mntd(metrics.input$picante.cdm, metrics.input$dists, abundance.weighted=FALSE)

aw_mntd <- function(metrics.input)
	mntd(metrics.input$picante.cdm, metrics.input$dists, abundance.weighted=TRUE)

my_psv <- function(metrics.input)
{
	PSV <- psv(metrics.input$picante.cdm, metrics.input$tree)
	PSV <- PSV$PSVs
	PSV
}

my_psc <- function(metrics.input)
{
	PSC <- pscCorr(metrics.input$picante.cdm, metrics.input$tree)
	PSC <- PSC$PSCs
	PSC
}

my_pse <- function(metrics.input)
{
	PSE <- pse(metrics.input$picante.cdm, metrics.input$tree)
	PSE <- PSE$PSEs
	PSE
}

my_PD <- function(metrics.input)
{
	PD <- pd(metrics.input$picante.cdm, metrics.input$tree, include.root=TRUE)
	PD <- PD$PD
	PD
}

my_PD_Cadotte <- function(metrics.input)
{
	PD <- pd(metrics.input$picante.cdm, metrics.input$tree, include.root=FALSE)
	PD <- PD$PD
	PD
}

my_QE <- function(metrics.input)
{
	QE <- raoD(metrics.input$picante.cdm, metrics.input$tree)
	QE <- QE$Dkk
	QE
}

##########################################################################################
################################### BETA METRICS #########################################
##########################################################################################

my_betaRichness <- function(metrics.input)
	sum(apply(metrics.input$picante.cdm, 2, lengthNonZeros) != 0)

my_totalAbundance <- function(metrics.input)
	sum(metrics.input$picante.cdm)

#should consider adding a skewness function here. need to check, but I believe the formula
#may just be (3(mean-median))/standard deviation

my_Ist <- function(metrics.input)
	spacodi.calc(sp.plot=t(metrics.input$picante.cdm), phy=metrics.input$tree)$Ist

my_Pst <- function(metrics.input)
	spacodi.calc(sp.plot=t(metrics.input$picante.cdm), phy=metrics.input$tree)$Pst

my_Bst <- function(metrics.input)
	spacodi.calc(sp.plot=t(metrics.input$picante.cdm), phy=metrics.input$tree)$Bst

my_PIst <- function(metrics.input)
	spacodi.calc(sp.plot=t(metrics.input$picante.cdm), phy=metrics.input$tree)$PIst

mean_naw_mpd <- function(metrics.input)
	mean(modifiedMPD(metrics.input$picante.cdm, metrics.input$dists,
	abundance.weighted=FALSE))

mean_inter_mpd <- function(metrics.input)
	mean(modifiedMPD(metrics.input$picante.cdm, metrics.input$dists,
	abundance.weighted="interspecific"))

mean_intra_mpd <- function(metrics.input)
	mean(modifiedMPD(metrics.input$picante.cdm, metrics.input$dists,
	abundance.weighted="intraspecific"))

mean_complete_mpd <- function(metrics.input)
	mean(modifiedMPD(metrics.input$picante.cdm, metrics.input$dists,
	abundance.weighted="complete"))

mean_naw_mntd <- function(metrics.input)
	mean(mntd(metrics.input$picante.cdm, metrics.input$dists,
	abundance.weighted=FALSE))

mean_aw_mntd <- function(metrics.input)
	mean(mntd(metrics.input$picante.cdm, metrics.input$dists,
	abundance.weighted=TRUE))

##########################################################################################
############################### ERROR CHECKING WRAPPERS ##################################
##########################################################################################

# lapply wrapper for Wilcoxon-signed rank tests
#
# Just a utility function. lapplies wilco.test functions over a list of dataframes.
#
# #param null.list A list of dataframes, one per null model, of observed metric scores.
# #param alternative Optional alternative hypothesis. Default is "two-sided". Use 
# "greater" for competition; "less" for habitat filtering.
#
# #details lapplies wilcoWrapApply over a list of dataframes. 
#
# #return A dataframe, with one row for each metric. The first column is the mean of the
# vector of metric values, the second is the p.value of whether it differs from mu=0,
# and the third is the name of the metric.
#
# #references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
# structure metrics and null models: a review with new methods and software.
# bioRxiv 025726.
#
# #examples
# a <- rnorm(n=100)
# b <- rnorm(n=100, mean=100)
# ex <- data.frame(a, b)
# test <- list("ex1"=ex, "ex2"=ex)
# wilcoWrapLApply(test, alternative="two.sided")

wilcoWrapLApply <- function(null.list, alternative)
{
	#lapply wilcoWrapApply over null.list. note that null.list is a list of data frames
	#from a single spatial simulation. you will further lapply this function over a list
	#of lists of data frames in a final test function
	temp <- lapply(null.list, wilcoWrapApply, alternative=alternative)

	#reduce the output list into a single data frame
	output <- Reduce(rbind, temp)

	#create a vector of expanded null model names. note that this code is sensitive to
	#changes. for instance, if one null model tests certain metrics that another does not,
	#this will not end up being correct. this generates a data frame, but we only save the
	#first column
	nullNames <- expand.grid(temp[[1]]$metric, names(null.list))[,2]
	
	output$null.model <- nullNames
	
	output
}

wilcoWrapApply <- function(dataframe, alternative)
{
	#exclude "richness" and "plot" columns
	exclude <- c("richness", "plot")
	temp <- dataframe[ ,!(names(dataframe) %in% exclude)]

	#apply wilcoWrap over data frame of metric SES scores for a given null and spatial sim
	output <- apply(temp, 2, wilcoWrap, mu=0, alternative)
	
	#transform the table, convert to a data frame, save the row names as an actual column,
	#exclude "richness" as a metric. output a data frame with three columns
	output <- t(output)

	#convert to data frame
	output <- as.data.frame(output)
	
	#add column names
	names(output) <- c("estimate", "p.value")
	
	#get rid of row names
	output$metric <- row.names(output)
	
	row.names(output) <- NULL

	output
}

wilcoWrap <- function(vect, mu=0, alternative)
{
	if(missing(alternative))
	{
		alternative <- "two.sided"
	}

	#set up a blank matrix to save results into
	output <- matrix(nrow=1, ncol=2)
	
	#run a quick t.test on the vector
	temp <- wilcox.test(x=vect, mu=mu, alternative=alternative)

	#pull out the observed mean and p.value from temp and retain these
	output[1,1] <- mean(vect, na.rm=TRUE)
	output[1,2] <- temp$p.value
	
	output
}

# lapply wrapper for t-tests
#
# Just a utility function. lapplies t-tests over lists of data frames
#
# #param null.list A list of dataframes, one per null model, of observed metric scores.
#
# #details lapplies tWrapApply over a list of dataframes. There is currently no easy
# way to pass arguments (like changing mu or whether the test is two- or one-tailed) down
# to tWrap.
#
# #return A dataframe, with one row for each metric. The first column is the mean of the
# vector of metric values, the second is the p.value of whether it differs from mu=0,
# and the third is the name of the metric.
#
# #export
#
# #references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
# structure metrics and null models: a review with new methods and software.
# bioRxiv 025726.
#
# #examples
# a <- rnorm(n=100)
# b <- rnorm(n=100, mean=100)
# ex <- data.frame(a, b)
# test <- list("ex1"=ex, "ex2"=ex)
# tWrapLApply(test)

tWrapLApply <- function(null.list)
{
	#lapply tWrapApply over null.list
	temp <- lapply(null.list, tWrapApply)

	#reduce the output list into a single data frame
	output <- Reduce(rbind, temp)

	#create a vector of expanded null model names. note that this code is sensitive to
	#changes. for instance, if one null model tests certain metrics that another does not,
	#this will not end up being correct. this generates a data frame, but we only save the
	#first column
	nullNames <- expand.grid(temp[[1]]$metric, names(null.list))[,2]
	
	output$null.model <- nullNames
	
	output
}

tWrapApply <- function(dataframe)
{
	#exclude "richness" and "plot" columns
	exclude <- c("richness", "plot")
	temp <- dataframe[ ,!(names(dataframe) %in% exclude)]

	#apply tWrap over a data frame of metric SES scores for a given null and spatial sim
	output <- apply(temp, 2, tWrap)
	
	#transform the table, convert to a data frame, save the row names as an actual column,
	#exclude "richness" as a metric. output a data frame with three columns
	output <- t(output)

	#convert to data frame
	output <- as.data.frame(output)
	
	#add column names
	names(output) <- c("estimate", "p.value")
	
	#get rid of row names
	output$metric <- row.names(output)
	
	row.names(output) <- NULL

	output
}

tWrap <- function(vect, mu=0)
{
	#set up a blank matrix to save results into
	output <- matrix(nrow=1, ncol=2)
	
	#run a quick t.test on the vector
	temp <- t.test(x=vect, mu=mu)
	
	#pull out the observed mean and p.value from temp and retain these
	output[1,1] <- temp$estimate
	output[1,2] <- temp$p.value
	
	output
}

##########################################################################################
################################ OTHER MISC FUNCTIONS ####################################
##########################################################################################

# Utility function to identify minimum values
#
# Given a vector where the last element is the minimum, identifies which elements in that
# vector match the last element.
#
# #param x A vector
# 
# #details Simple utility function, used in the competitive exclusion simulations.
#
# #return A logical vector of length input vector minus 1, corresponding to whether an
# element of the input vector equals the last element of the input vector.
#
# #references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
# structure metrics and null models: a review with new methods and software.
# bioRxiv 025726.
#
# #examples
# #create a basic input vector
# temp <- c(1,2,3,4,5,6,1)
#
# #here is an example of the compareMins function 
# compareMins(temp)

compareMins <- function(x)
{
	output <- x[1:(length(x)-1)] == x[length(x)]
	return(output)
}

# Scale output of evolveTraits to arena size
#
# Given a matrix of two traits, and the minimum and maximum extent of the desired arena,
# will return a data frame of species' traits scaled to the arena size.
#
# #param input.traits Second element of the results of a call to evolveTraits()
# #param min.arena Minimum size of arena, e.g. 0
# #param max.arena Maximum size of arena
# 
# #details Scales a matrix of species' traits to a desired mininimum-maximum range.
# Intended for use in a spatially explicit scenario with two traits, but could easily
# be co-opted.
#
# #return A scaled and named dataframe of species traits
#
# #export
#
# #references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
# structure metrics and null models: a review with new methods and software.
# bioRxiv 025726.
#
# #examples
# tree <- geiger::sim.bdtree(b=0.1, d=0, stop="taxa", n=50)
#
# temp <- evolveTraits(tree)
#
# scaled <- scaler(temp[[2]], min.arena=0, max.arena=300)

##write a function that will take the second element of the output of the evolveTraits
##function, and the min and max arena arguments, and output a data frame of scaled traits
##where min and max traits are min and max of arena

scaler <- function(input.traits, min.arena, max.arena)
{
	std1 <- (input.traits[,1] - 
		min(input.traits[,1]))/(max(input.traits[,1])-min(input.traits[,1]))
	std2 <- (input.traits[,2] - 
		min(input.traits[,2]))/(max(input.traits[,2])-min(input.traits[,2]))
	
	output.trait1 <- (max.arena - min.arena) * std1 + min.arena
	output.trait2 <- (max.arena - min.arena) * std2 + min.arena
	
	output.traits <- cbind(output.trait1, output.trait2)
	
	return(output.traits)
}

# Identify failed runs
#
# This identifies spatial sim/null model runs that failed
#
# #param single.iteration The results of a single iteration of multiLinker
# #param concat.by Whether randomizations were concatenated by richness, plot or both.
#
# #details It is possible that a given null model, e.g. regional, failed on a given run. 
# These failures, in our experience so far, can only occur if a given species richness is
# insufficiently sampled by the null model when concat.by="richness". Otherwise, I
# believe the run will always succeed. Thus, to simplify coding, this function only
# evaluates the by "richness" runs for failure when concat.by="both". If it is possible
# for a null to fail if concatenating by plot then this function will need to be
# re-written. The main point of this function it to identify failed runs and later
# remove them from the results of a given iteration. This prevents errors further
# down the line when summarizing results.
#
# #return A data frame summarizing which simulation/null runs failed.
#
# #export
#
# #references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
# structure metrics and null models: a review with new methods and software.
# bioRxiv 025726.
#
# #examples
# #not run
# #results <- readIn()
# #summ <- failed(results[[1]], "both")

failed <- function(single.iteration, concat.by)
{
	#it is really tough to properly go across each sub element within a single iteration
	#result with lapply, so just get in there with some for loops and brute force it.
	#First set up an empty matrix where we will save the dimensions of each null model
	#from each spatial simulation. If any of these have 0s, then something failed with
	#the null model. Similarly, the regional null very infrequently fails and throws NAs.
	#For these, if the sum of is.na of any the columns equals the length of the column
	#then there is a problem. So,
	#add a part at end of the for loop where even if the dimensions are ok, if the sum of
	#NAs equals the length then we say the dimensions are 0. This will allow easy 
	#identification later. The following dimensions assume that each spatial simulation 
	#was tested with the same null models and each has an element named ses.
	
	temp <- matrix(ncol=4, nrow=length(single.iteration) * 
		length(single.iteration[[1]]$ses))
	
	temp <- as.data.frame(temp)
	
	for(i in 1:length(single.iteration))
	{
		for(j in 1:length(single.iteration[[i]]$ses))
		{
			#define the length of the null models
			nullLength <- length(single.iteration[[i]]$ses)
			#now define the row ID as nullLength times i minus 1 + j
			rowID <- nullLength * (i - 1) + j
			#pull the relevant name for the spatial simulation
			temp[rowID,1] <- names(single.iteration)[i]
			#pull the relevant name for the null model
			temp[rowID,2] <- names(single.iteration[[i]]$ses)[j]
			#find the dimensions of the null model
			if(concat.by=="richness" | concat.by=="plot")
			{
				dims <- dim(single.iteration[[i]]$ses[[j]])
			}
			else if(concat.by=="both")
			{
				dims <- dim(single.iteration[[i]]$ses[[j]]$richness)
			}
			else
			{
				stop("concat.by must equal either both, richness, or plot")
			}
			#set the relevant row equal to the dimensions
			temp[rowID,3:4] <- dims
			#now add a cheap fix to skip to next iteration of loop if the dimensions have
			#0 or 1 or else the lengthNA double apply statement below will fail
			if(temp[rowID,3] <= 1)
			{
				next()
			}
			#otherwise go into the remainder of the for loop
			#determine the length of each column and of NAs per column
			if(concat.by=="richness" | concat.by=="plot")
			{
				lengthNA <- apply(apply(single.iteration[[i]]$ses[[j]], 2, is.na), 2, sum)
				lengthCol <- apply(single.iteration[[i]]$ses[[j]], 2, length)
			}
			else
			{
				lengthNA <- apply(apply(single.iteration[[i]]$ses[[j]]$richness, 
					2, is.na), 2, sum)
				lengthCol <- apply(single.iteration[[i]]$ses[[j]]$richness, 2, length)
			}
			#delete the element named "richness" from each of these
			lengthNA <- lengthNA[names(lengthNA) != "richness"]
			lengthCol <- lengthCol[names(lengthCol) != "richness"]
			#now if these things are ever equal there is a problem
			if(sum(lengthNA == lengthCol) > 0)
			{
				temp[rowID,3:4] <- 0
			}
		}
	}
	
	#set the column names to be nice for output
	names(temp) <- c("simulation", "null", "dim1", "dim2")

	#subset it to instances where dim1 == 0.
	temp <- temp[temp$dim1 == 0,]
	
	temp
	
}

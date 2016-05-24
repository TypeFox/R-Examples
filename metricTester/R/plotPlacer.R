#' Randomly place plots in arena
#'
#' Given a desired number of plots, the arena size, and the plot size, will attempt
#' to place plots down in a non-overlapping fashion
#'
#' @param no.plots Number of plots to place
#' @param arena.length Length of one side of arena
#' @param plot.length Length of one side of desired plot. 
#' 
#' @details Places plots down in non-overlapping fashion according to parameters
#' supplied. Because this would run indefinitely if unacceptable parameters were supplied,
#' a conservative check is implemented to "ensure" the function does not get stuck. If
#' unacceptable parameters are supplied, will return an arena and a smaller total sampling
#' area will need to be defined.
#'
#' @return A matrix with the X & Y coordinates of the four corners of each plot placed
#'
#' @export
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#'
#' boundResults <- plotPlacer(no.plots=10, arena.length=300,
#'	plot.length=50)

plotPlacer <- function(no.plots, arena.length, plot.length)
{
	#this ugly bit of code is because I used to define these things as arguments to
	#function. this makes input easier. should fix code to just use arena.length and not
	#require x.max and y.max
	x.max=arena.length
	y.max=arena.length

	if(((plot.length^2)/(arena.length^2))*no.plots > 0.4)
	{
		stop("Plot and/or arena parameters unsuitable. Sample less of total arena")
	}

	##define plot bounds, etc

	plot.bounds <- matrix(0,nrow=no.plots,ncol=4)
	colnames(plot.bounds) <- c("X1","X2","Y1","Y2")

	for (i in c(1:no.plots))
	{
		repeat 
		{
			OK <- TRUE
			#sample a point that is within the uniform distribution from 0 to the arena
			#bounds minus space for the plot bounds (vector is 1000 long for now).
			plot.bounds[i,1] <- sample(runif(n=1000, min=0, 
				max=x.max-plot.length), 1)
			#add the plot length to this number. note that this means that X2 is always
			#more than X1
			plot.bounds[i,2] <- plot.bounds[i,1] + plot.length
			#as above
			plot.bounds[i,3] <- sample(runif(n=1000, min=0, 
				max=y.max-plot.length), 1)
			plot.bounds[i,4] <- plot.bounds[i,3] + plot.length
			if (i > 1) 
			{
				for (j in c(1:(i-1)))
				{
					#write an if statement where if any of the conditions are true, 
					#OK is set to FALSE and it loops back through
					if (any(
						#the corner X1,Y1 is within the bounds of another plot
						plot.bounds[i,1] > plot.bounds[j,1] & 
							plot.bounds[i,1] < plot.bounds[j,2] & 
							plot.bounds[i,3] > plot.bounds[j,3] &
							plot.bounds[i,3] < plot.bounds[j,4],
						#the corner X2,Y1 is within the bounds of another plot
						plot.bounds[i,2] > plot.bounds[j,1] & 
							plot.bounds[i,2] < plot.bounds[j,2] & 
							plot.bounds[i,3] > plot.bounds[j,3] &
							plot.bounds[i,3] < plot.bounds[j,4],
						#the corner X1,Y2 is within the bounds of another plot
						plot.bounds[i,1] > plot.bounds[j,1] & 
							plot.bounds[i,1] < plot.bounds[j,2] & 
							plot.bounds[i,4] > plot.bounds[j,3] &
							plot.bounds[i,4] < plot.bounds[j,4],
						#the corner X2,Y2 is within the bounds of another plot
						plot.bounds[i,2] > plot.bounds[j,1] & 
							plot.bounds[i,2] < plot.bounds[j,2] & 
							plot.bounds[i,4] > plot.bounds[j,3] &
							plot.bounds[i,4] < plot.bounds[j,4]
						)) 
					{
						OK <- FALSE
					}
				}
			}
			if (OK == TRUE) 
			{
				break;
			}
		}
	}
	x.centers <- apply(plot.bounds[,1:2], 1, mean)
	y.centers <- apply(plot.bounds[,3:4], 1, mean)
	
	centers <- data.frame(x.center=x.centers, y.center=y.centers)

	#give the centers names so you can later subset it to those that had sufficient
	#individuals in based on plotContents
	row.names(centers) <- paste("tempPlot", 1:no.plots, sep="")
	dists <- as.matrix(dist(centers, diag=TRUE, upper=TRUE))

	output <- list("dists"=dists, "plot.bounds"=plot.bounds)

	return(output)
}

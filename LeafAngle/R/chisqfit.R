`chisqfit` <-
function(angles, distribution, distpars=NA, ...){

		# set the breaks
		stepsize <- 5
		lbreaks <- seq(0,90,by=stepsize)
		midpoints <- lbreaks + stepsize/2
		midpoints <- midpoints[-length(midpoints)]
		N <- length(angles)
		
		# Expected counts based on true distribution.
		# Density function is evaluated in the middle of the bin
		# (Should really be integrated, but this is good enough).
		fdens <- ftheta(midpoints, degrees=TRUE, distribution=distribution, distpars=distpars, ...)
		binprobs <- stepsize * fdens
		simcounts <- binprobs * N
		
		# Histogram of the data
		hdata <- hist(angles, breaks=lbreaks, plot=FALSE)
		datacounts <- hdata$counts	
		
	    # Take only the counts where they are > 0.
		part1 <- (datacounts - simcounts)^2
		part1 <- part1[simcounts > 0]
		simcounts <- simcounts[simcounts > 0]
			
		chisq <- sum(part1  / simcounts ) 

return(chisq)
}


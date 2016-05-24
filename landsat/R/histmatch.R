histmatch <-
function(master, tofix, mask, minval=0, maxval=255, by=1)
{
	# simple histogram matching function
    # mask should contain NA for values to use; all other values will be omitted
	results <- tofix # want to return results in same format
	master <- as.vector(as.matrix(master))
	tofix <- as.vector(as.matrix(tofix))

    if(missing(mask)) mask <- rep(NA, length(master))
    else mask <- as.vector(as.matrix(mask))
    results.final <- rep(NA, length(mask))

    master <- master[is.na(mask)]
    tofix <- tofix[is.na(mask)]

	breaks <- seq(minval, maxval, by=by)
	master.cdf <- hist(master, breaks=breaks, plot=FALSE) 
	master.cdf <- c(0, cumsum(master.cdf$counts/sum(master.cdf$counts)))
	tofix.cdf <- hist(tofix, breaks=breaks, plot=FALSE) 
	tofix.cdf <- c(0, cumsum(tofix.cdf$counts/sum(tofix.cdf$counts)))

    # fixed 2012-07-16 to work with continuous data
    # originally written to work with integer data
	results.recode <- breaks
    results.values <- rep(NA, length(tofix))
    # original #	for(i in 1:length(breaks)) {
    # original #        testvals <- breaks[master.cdf < tofix.cdf[i]]
    # original #        if(length(testvals) > 0)
    # original #            results.recode[i] <- max(testvals)
    # original #        results.values[tofix == breaks[i]] <- results.recode[i]
    # original #    }

    for (i in 2:length(breaks)) {
        testvals <- breaks[master.cdf < tofix.cdf[i]]
        if (length(testvals) > 0) 
            results.recode[i] <- max(testvals)
        results.values[tofix > breaks[i-1] & tofix <= breaks[i]] <- results.recode[i]
    }

    results.final[is.na(mask)] <- results.values

    if(class(results) == "SpatialGridDataFrame")
        results@data[,1] <- results.final
    else if(is.data.frame(results))
        results <- data.frame(matrix(results.final, nrow=nrow(results), ncol=ncol(results)))
    else if(is.matrix(results))
        results <- matrix(results.final, nrow=nrow(results), ncol=ncol(results))
    else
        results <- results.final

    list(recode=results.recode, newimage=results)
}


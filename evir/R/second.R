"decluster" <- 
function(series, run = NA, picture = TRUE)
{
    n <- length(as.numeric(series))
    times <- attributes(series)$times
    if(is.null(times)) stop("`series' must have a `times' attribute")
    as.posix <- is.character(times) || inherits(times, "POSIXt") ||
      inherits(times, "date") || inherits(times, "dates")
    if(as.posix)
        gaps <- as.numeric(difftime(as.POSIXlt(times)[2:n],
            as.POSIXlt(times)[1:(n-1)], units = "days"))
    else gaps <- as.numeric(diff(times))
    longgaps <- gaps > run
    if(sum(longgaps) <= 1)
        stop("Decluster parameter too large")
    cluster <- c(0, cumsum(longgaps))
    whichcmax <- tapply(as.numeric(series), cluster, which.max)
    clen <- cumsum(tapply(as.numeric(series), cluster, length))
    whichcmax <- whichcmax + c(0, clen[-length(clen)])
    cmax <- as.numeric(series)[whichcmax]
    newtimes <- times[whichcmax]
    newseries <- structure(series[whichcmax], times = newtimes)
    n <- length(as.numeric(newseries))

    if(as.posix) {
        newgaps <- as.numeric(difftime(as.POSIXlt(newtimes)[2:n],
            as.POSIXlt(newtimes)[1:(n-1)], units = "days"))
        times <- as.POSIXlt(times)
        newtimes <- as.POSIXlt(newtimes)
    }
    else newgaps <- as.numeric(diff(newtimes))
    
    if(picture) {
      	cat("Declustering picture...\n")
       	cat(paste("Data reduced from", length(as.numeric(series)),
       		"to", length(as.numeric(newseries)), "\n"))
       	par(mfrow = c(2, 2))
        plot(times, series, type = "h")
	qplot(gaps)
        plot(newtimes, newseries, type = "h")
       	qplot(newgaps)
       	par(mfrow = c(1, 1))
    }
    newseries
}

"findthresh" <- 
function(data, ne)
{
	data <- rev(sort(as.numeric(data)))
	thresholds <- unique(data)
	indices <- match(data[ne], thresholds)
	indices <- pmin(indices + 1, length(thresholds))
	thresholds[indices]
}


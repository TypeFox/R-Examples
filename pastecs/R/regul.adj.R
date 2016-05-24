"regul.adj" <-
function(x, xmin=min(x), frequency=NULL, deltat=(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))/(length(x)-1), tol=deltat, tol.type="both", nclass=50, col=c(4, 5, 2), xlab=paste("Time distance"), ylab=paste("Frequency"), main="Number of matching observations", plotit=TRUE, ...) {
	xmin <- xmin
	frequency <- frequency
	deltat <- deltat
	x <- sort(x)				# Make sure x is sorted in increasing order
	x <- x[!is.na(x)]			# Eliminate missing values
	n <- length(x)
	if (is.null(frequency)) {
		if (is.null(deltat)) {
			stop("You must define at least one of frequency or deltat")
		} else {
			frequency <- 1/deltat
		}
	} else {
		deltat <- 1/frequency
	}
	# We verify also that tol is a round fraction of deltat
	if (is.null(tol) || tol == 0) tol2 <- 0 else {
		tol2 <- abs(tol)
		if (tol2 > deltat) tol2 <- deltat else {
			tol2 <- deltat/round(deltat/tol2)
		}
	}
	# We calculate highest n, so as there is no extrapolation
	if (max(x) < xmin) {		# There is no solution!
		nx <- 0
		stop("xmin is higher than max value in x!")
	} else {					# We calculate nx
		nx <- floor((max(x) - xmin)/deltat) + 1
	}
	# create xout vector
	xout <- 0:(nx-1) * deltat + xmin
	# calculate the matching vector
	pos <- match.tol(xout, x, tol.type=tol.type, tol=tol2)
	# which regular date match observations? Put distance for a matching observation, Inf for an interpolated observation, -Inf for an extrapolated observation
	match.dist <- xout - x[pos]
	match.dist[is.na(match.dist)] <- Inf
	match.dist[(xout < min(x) | xout > max(x)) & match.dist == Inf] <- -Inf
	# Number of matches
	match <- sum(is.finite(match.dist))
	exact.match <- sum(match.dist == 0)
	# Construct the params vector
	params <- c(xmin, nx, deltat, tol2)
	names(params) <- c("xmin", "n", "deltat", "tol")
	# Draw the graph, if plot is TRUE
	if (plotit == TRUE) {
		if (tol2 == 0) HT <- 1.001 else HT <- 101*tol2/100
		Data <- abs(match.dist)
		Data[is.infinite(Data)] <- HT 			# Inf are replaced by a value higher than Tol
		Data[Data == 0] <- -0.00001				# For putting exact matching values in a separate category
		# Don't draw, but get vectors of results
		res <- hist(Data, nclass=nclass, plot=FALSE)
		classes <- res$breaks[2:length(res$breaks)]
		ncl <- length(classes)
		classes[ncl] <- Inf
		counts <- res$counts
		names(counts) <- classes
		# Create a vector for colors, so as the first and last classes are drawn in a different color
		cols <- NULL
		cols[1] <- col[2]
		if (sum(Data == -0.00001) > 0) cols[1] <- col[1]
		if (ncl > 2) cols[2:(ncl-1)] <- col[2]
		cols[ncl] <- col[3]
		# Actually draw the histogram
		hist(Data, nclass=nclass, col=cols, xlab=xlab, ylab=ylab, main=main)
		counts <- counts[counts != 0]
		lc <- length(counts)
		counts2 <- NULL
		for (i in 1:lc) {
			counts2[i] <- sum(counts[1:i])
		}
		names(counts2) <- names(counts)
		res <- list(params=params, match=match, exact.match=exact.match, match.counts=counts2)
	} else {
		res <- list(params=params, match=match, exact.match=exact.match)
	}
	res
}

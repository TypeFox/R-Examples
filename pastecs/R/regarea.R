"regarea" <-
function (x, y = NULL, xmin=min(x), n=length(x), deltat=(max(x)-min(x))/(n-1), rule=1, window=deltat, interp=FALSE, split=100) {
	x <- xy.coords(x, y)
	y <- x$y
	x <- x$x
	# Make sure data are sorted in increasing order according to x
	srt <- sort.list(x)
	x <- x[srt]
	y <- y[srt]
	# Check arguments
	if (!is.numeric(x) || !is.numeric(y)) 
		stop("regarea: x and y must be numeric")
    nx <- length(x)
    if (nx != length(y)) 
        stop("x and y must have equal lengths")
    if (nx < 2) 
        stop("regarea requires at least two values to interpolate")
    # Eliminate entries with missing values
    ok <- !(is.na(x) | is.na(y))
    x <- x[ok]
    y <- y[ok]
    nx <- length(x)
    if (nx < 2) 
        stop("regarea requires at least two non-missing values to interpolate")
    if (window <= 0)
	    stop("the window must be a positive number")
	# The next tree lines must be activated if one want at least one observation INSIDE each window (but not necessary)
	#largestgap <- max(x[2:nx]-x[1:(nx-1)])
	#if (window <= largestgap)		# The window must contain at least one value everywhere!
    #	stop(paste("the window must be wider than the largest gap in the series (", largestgap, ")", sep=""))
    if (n <= 0) 
    	stop("regarea requires n >= 1")
    if (deltat <= 0) 
    	stop("regarea requires deltat > 0")
    xout <- 0:(n-1) * deltat + xmin							# vector of xout regular sequence
    xmax <- xout[n]
    if (!is.numeric(split)) 
		stop("regarea: split must be numeric")
	split <- round(split)
	if (split < 1) split <- 1								# split must be a positive integer!
    # Misc values required for calculations
    halfwin <- window/2
    xlb <- xmin - halfwin - 1								# first area must extend at least below lowest window
    xlb <- min(xlb, x[1])
    xub <- xmax + halfwin + 1								# last area must extend at least above highest window
    xub <- max(xub, x[nx])
	# all x must be positive => shift the scale to meet this requirement
    if (xlb <= 0) {
    	shift <- (-xlb + 1)
    	xout <- xout + shift
    	x <- x + shift
    	xlb <- xlb + shift
    	xub <- xub + shift
    } else shift <- 0
    # vector containing lower bound (x1i) of areas for each data xi
    x1 <- NULL
    x1[1] <- xlb											# Make sure the first area include (xmin - window/2)
    x1[2:nx] <- (x[1:(nx-1)]+x[2:nx])/2						# Areas start in the middle of two data
    # vector containing upper bound (x2i) of areas for each data xi
    x2 <- NULL
    x2[1:(nx-1)] <- x1[2:nx]								# One upper bound is equal to the next lower bond
    x2[nx] <- xub											# Make sure the last area include (xmax + window/2)
    # matrix containing lower bound of each window for xout
    u1 <- xout - halfwin
	# To avoid too large matrices, and a lot of wasted calculations, the series is split is shorter sub-series
    # that are analyzed sequencially
    seq <- ceiling(n/split)
    yout <- NULL
    for (i in 1:seq) {
    	# we select submatrices corresponding to the current sequence to analyse
    	sn <- split
    	sl <- (i - 1) * split + 1
    	su <- i * split
    	if (su > n) {										# When the last sequence is partial only (less than n elements remains)
    		su <- n
    		sn <- su - sl + 1
    	}
    	smin <- xout[sl] - halfwin
    	smax <- xout[su] + halfwin
    	slx <- sum(x1 < smin)
    	sux <- sum(x2 < smax) + 1
    	snx <- sux - slx + 1
    	# submatrices xs1 and xs2 are extracted and duplicated sn times accross columns
    	xs1 <- rep(x1[slx:sux], sn)
     	dim(xs1) <- c(snx, sn)
    	xs2 <- rep(x2[slx:sux], sn)
		dim(xs2) <- c(snx, sn)
    	# submatrices us1 and us2 are extracted and duplicated across snx rows
    	us1 <- rep(u1[sl:su], snx)
    	dim(us1) <- c(sn, snx)
    	us1 <- t(us1)
    	us2 <- us1 + window
    	# keep only maximal value of xs1 or us1
    	us1 <- pmax(xs1, us1)
    	# keep only minimal value of xs2 or us2
    	us2 <- pmin(xs2, us2)
    	# fi = u2i - u1i
    	xs1 <- us2 - us1
    	# negative values for f indicate areas totally outside the window => f is reset to 0 for them
    	xs1[xs1<0] <- 0
    	# approximation of yout (for the submatrix) is the sum of areas (y * f, in columns) divided by the window length
    	yout[sl:su] <- apply((y[slx:sux] * xs1), 2, sum)/window
    }
    # If interp is false, we replace matching values of yout by y for each xout that match an x
    if (interp == FALSE) {
    	pos <- match(xout, x, nomatch=NA)
    	youtmatch <- y[pos]											# Any matching value in the right place in yout, the rest is NAs
    	posmatch <- !is.na(pos)
    	yout[posmatch] <- youtmatch[posmatch]						# and they replace extrapolated values in yout
    }
    # If rule == 1, we still must replace values outside {x[1], x[nx]} by NA
    if (rule == 1)
    	yout[xout < x[1] | xout > x[nx]] <- NA
    # Restore initial xout values
    xout <- xout - shift
    # And finally, we return the list of reguled data (x, y)
    list(x = xout, y = yout)
}

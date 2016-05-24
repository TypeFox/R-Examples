plot.hwtANYN <-
function (x, xlabvals, xlabchars, ylabchars, first.level = 1, 
    main = "Haar Wavelet Coefficients", scaling = c("global", "by.level"), 
    rhlab = FALSE, sub, NotPlotVal = 0.005, xlab = "Translate", 
    ylab = "wd-equivalent Resolution Level", miss.coef.col=2, miss.coef.cex=0.5,
    miss.coef.pch=2, ...) 
{
    scaling <- match.arg(scaling)
    ctmp <- class(x)
    if (is.null(ctmp)) 
        stop("hwtANYN has no class")
    else if (ctmp != "hwtANYN") 
        stop("hwtANYN is not of class hwtANYN")

    if (x$reindex==FALSE)
	stop("Can only plot object with reindex=TRUE")

    if (first.level < 1)
	stop("first.level has to be >= 1")

    levels <- nlevelsWT(x)

    if (first.level >= levels)
	stop(paste("first.level has to be < ", levels))


    nlevels <- levels - first.level

    type <- x$type

    if (type == "wavelet") 
        n <- 2^(levels - 1)
    else if (type == "station") 
        n <- 2^levels
    else stop("Unknown type for wavelet object")

    if (missing(sub)) 
        sub <- paste(switch(type, wavelet = "Standard transform", 
            station = "Nondecimated transform"),"Haar")

    plot(c(0, 0, n, n), c(0, nlevels + 1, nlevels + 1, 0), type = "n", 
        xlab = xlab, ylab = ylab, main = main, yaxt = "n", xaxt = "n", 
        sub = sub, ...)

    yll <- (levels - 1):first.level

    if (missing(ylabchars)) 
        axis(2, at = 1:(nlevels), labels = yll)
    else if (length(ylabchars) != nlevels) 
        stop(paste("Should have ", nlevels, " entries in ylabchars"))
    else axis(2, at = 1:(nlevels), labels = ylabchars)

    if (missing(xlabchars)) {
        if (missing(xlabvals)) {
            if (type == "wavelet") 
                axx <- c(0, 2^(levels - 3), 2^(levels - 2), 2^(levels - 
                  2) + 2^(levels - 3), 2^(levels - 1))
            else axx <- c(0, 2^(levels - 2), 2^(levels - 1), 
                2^(levels - 1) + 2^(levels - 2), 2^levels)
            if (is.null(tsp(x))) 
                axis(1, at = axx)
            else {
                v <- seq(from = tsp(x)["start"], by = tsp(x)["deltat"], 
                  length = n)
                if (type == "wavelet") 
                  atl <- 2 * v
                else atl <- v
                atl <- pretty(atl, n = 4)
                ats <- (n * atl)/(max(atl) - min(atl))
                axis(1, at = ats, labels = atl)
            }
        }
        else {
            lx <- pretty(xlabvals, n = 4)
            cat("lx is ", lx, "\n")
            if (lx[1] < min(xlabvals)) 
                lx[1] <- min(xlabvals)
            if (lx[length(lx)] > max(xlabvals)) 
                lx[length(lx)] <- max(xlabvals)
            cat("lx is ", lx, "\n")
            xix <- NULL
            for (i in 1:length(lx)) {
                u <- (xlabvals - lx[i])^2
                xix <- c(xix, (1:length(u))[u == min(u)])
            }
            axx <- xix
            if (type == "wavelet") 
                axx <- xix/2
            axl <- signif(lx, digits = 2)
            axis(1, at = axx, labels = axl)
        }
    }
    else axis(1, at = xlabvals, labels = xlabchars)
    myxx <- 1:n
    height <- 1
    axr <- NULL

    if (scaling == "global") {
        my <- 0
        for (i in (levels:first.level)) {
            y <- x$d[[i]]
	    if (!is.null(y))	{
		    y <- y[!is.na(y)]
		    my <- max(c(my, abs(y)))
		    }
        }
	cat("my is: ", my, "\n")
    }

    for (i in (levels:(first.level+1))) {
	y <- x$d[[i]]

	if (!is.null(y))	{

		y <- y[!is.na(y)]

		if (type == "wavelet") 
		    n <- 2^(i-1)
		else {
		    n <- 2^levels
		}
		xplot <- myxx
		ly <- length(y)
		if (scaling == "by.level") 
		    my <- max(abs(y))
		if (my == 0) {
		    y <- rep(0, length(y))
		}
		else y <- (0.5 * y)/my
		axr <- c(axr, my)
		if (max(abs(y)) > NotPlotVal) {
		    segments(1:ly, height, 1:ly, height + y)
		    }
		if (ly < n)
			points((ly+1):n, rep(height, n-ly), pch=miss.coef.pch,
				col=miss.coef.col, cex=miss.coef.cex)

		if (i != first.level) {
		    if (type == "wavelet") {
			x1 <- myxx[seq(1, n - 1, 2)]
			x2 <- myxx[seq(2, n, 2)]
			myxx <- (x1 + x2)/2
		    }
		    height <- height + 1
		}
	}
    }
    if (rhlab == TRUE) 
        axis(4, at = 1:length(axr), labels = signif(axr, digits = 3))
    axr
}

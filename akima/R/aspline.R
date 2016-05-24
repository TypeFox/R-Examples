aspline <- function(x, y=NULL, xout, n = 50, ties=mean, method="original", degree=3) {

    if (! method %in% c("original", "improved")) stop(paste("unknown method:", method))

    x <- xy.coords(x, y) # -> (x,y) numeric of same length
    y <- x$y
    x <- x$x
    nx <- length(x)
    if(any(na <- is.na(x) | is.na(y))) {
	ok <- !na
	x <- x[ok]
	y <- y[ok]
	nx <- length(x)
    }
    if (!identical(ties, "ordered")) {
	if (length(ux <- unique(x)) < nx) {
	    if (missing(ties))
		warning("Collapsing to unique x values")
	    y <- as.vector(tapply(y,x,ties))# as.v: drop dim & dimn.
	    x <- sort(ux)
	    nx <- length(x)
	} else {
	    o <- order(x)
	    x <- x[o]
	    y <- y[o]
	}
    }
    if (nx <= 1)
	stop("need at least two non-NA values to interpolate")
    if (missing(xout)) {
	if (n <= 0)
	    stop("aspline requires n >= 1")
	xout <- seq(x[1], x[nx], length = n)
    }
    nout <- length(xout)
    yout <- numeric(nout)
    err <- 0
    if (method == "improved") {
      ret <- .Fortran("uvip3p", as.integer(degree),
                    as.integer(nx), as.double(x), as.double(y), 
                    as.integer(nout), as.double(xout), yout=as.double(yout), 
                    err=as.integer(err), PACKAGE="akima")
    } else {
      ret <- .Fortran("intrpl",
                    as.integer(nx), as.double(x), as.double(y), 
                    as.integer(nout), as.double(xout), yout=as.double(yout), 
                    err=as.integer(err), PACKAGE="akima")
    }
    err <- max(0, min(ret$err, 10)) # 10 is maximum error code
    if (err > 0) {
      ## if the error handling befor .Fortran is correct
      ## the following lines should never be called
      errl <- c("Insufficient number of data points",
                "No desired points",
                "Two data points identical or out of sequence",
                "undefined error code",
                "undefined error code",
                "identical x values",
                "x values out of sequence",
                "undefined error code",
                "undefined error code",
                "internal error in Akima Fortran code")
      warning(errl[err])
    }
    list(x = xout, y = ret$yout)
}

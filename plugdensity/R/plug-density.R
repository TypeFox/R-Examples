## Idea : Change this to only estimate the bandwidth !
## ----> incorporate into base R

plugin.density <- function(x, nout = 201, xout = NULL, na.rm = FALSE)
{
    ## Purpose:	 Plug-in density estimate (global bandwidth)
    ## -------------------------------------------------------------------------
    ## Arguments: x: data;
    ##		nout : how many output values -> used if xout is NULL (default)
    ##		xout : explicit output abscissae
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 16 Mar 1998, 18:55

    if (!is.numeric(x))
	stop("argument must be numeric")
    name <- deparse(substitute(x))
    x <- as.vector(x)
    if (any(x.na <- is.na(x))) {
	if (na.rm)
	    x <- x[!x.na]
	else stop("x contains missing values")
    }
    n <- length(x <- sort(x))
    if(is.null(xout)) {
	## R's density() here extends the range depending on bandwidth !
	dx <- diff(rx <- range(x))
	if(dx < sqrt(.Machine$double.eps)) dx <- mean(abs(rx))/1000
	m <- as.integer(nout)
	xout <- seq(from=rx[1] - dx/10, to=rx[2] + dx/10, length= m)
    } else {
	m <- length(xout)
	if(is.unsorted(xout)) xout <- sort(xout)
    }
    r <- .C(plugin,
	    x = as.double(x), n=n,
	    z = xout, m=m,
	    f = double(m),
	    h = double(1))
    structure(list(x= r$z, y= r$f, bw = r$h, n=n,
		   call = match.call(), data.name = name),
	      class=c("densityEHpi", "density")) # Eva Herrman plug in
}

print.densityEHpi <- function(x, digits = getOption("digits"), ...)
{
    cat("EvaHerrmann plugin density estimate\n call :",
	deparse(x$call),"\n n = ", x$n,
	" ;  estimated (Gaussian) bandwidth h = ",
	format(x$bw, digits = digits),"\n")
    str(x[1:2], digits = digits, ...)
    invisible(x)
}

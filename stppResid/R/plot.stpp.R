plot.stpp <- function(x, ..., pch = 1, asp = 1)
{
	if(!is.stpp(x))
		stop("x must be an object of type stpp")
	plot(x$x, x$y, xlim = x$xcoord, ylim = x$ycoord, pch = pch, asp = asp, xlab = "x", ylab = "y", ...)
}
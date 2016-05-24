plot.thinresid <- function(x, ..., pch = 1, asp = 1)
{
	if(!is.thinresid(x))
		stop("x must be an object of type thinresid")
	plot(x[[3]]$x, x[[3]]$y, xlim = x[[1]]$xcoord, ylim = x[[1]]$ycoord, pch = pch, asp = asp, xlab = "x", ylab = "y", ...)
}
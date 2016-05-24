plot.superthin <- function(x, ..., pch1 = 1, pch2 = 3, asp = 1)
{
	if(!is.superthin(x))
		stop("x must be an object of type superthin")
	plot(x[[5]]$x, x[[5]]$y, xlim = x[[1]]$xcoord, ylim = x[[1]]$ycoord, pch = pch1, asp = asp, xlab = "x", ylab = "y", ...)
	points(x[[6]]$x, x[[6]]$y, pch = pch1)
	points(x[[4]]$x, x[[4]]$y, pch = pch2)	
}
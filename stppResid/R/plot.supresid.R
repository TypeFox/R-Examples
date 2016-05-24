plot.supresid <- function(x, ..., pch1 = 1, pch2 = 3, asp = 1)
{
	if(!is.supresid(x))
		stop("x must be an object of type supresid")
	plot(x[[1]]$x, x[[1]]$y, xlim = x[[1]]$xcoord, ylim = x[[1]]$ycoord, pch = pch1, asp = asp, xlab = "x", ylab = "y", ...)
	points(x[[4]]$x, x[[4]]$y, pch = pch2)
}
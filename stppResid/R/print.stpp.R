print.stpp <- function(x, ...)
{
	cat("Space-time point process\n")
	cat("Points:\n")
	cat("x: ")
	print(x$x)
	cat("y: ")
	print(x$y)
	cat("t: ")
	print(x$t)
	cat("Space-time window:\n")
	cat("x-range: ")
	print(x$xcoord)
	cat("y-range: ")
	print(x$ycoord)
	cat("t-range: ")
	print(x$tcoord)
}
print.stwin <- function(x, ...)
{
	cat("Space-time window\n")
	cat("x-range: ")
	print(x$xcoord)
	cat("y-range: ")
	print(x$ycoord)
	cat("t-range: ")
	print(x$tcoord)	
}
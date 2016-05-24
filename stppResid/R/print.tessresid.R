print.tessresid <- function(x, ...)
{
	print(x$X)
	cat("Tessellation residuals:\n")	
	print(x$residuals)
	cat("Tile list:\n")
	print(x$tile.list)
}
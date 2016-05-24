print.tessdev <- function(x, ...)
{
	print(x$X)
	cat("Tessellation deviance residuals:\n")	
	print(x$residuals)
	cat("Tile list:\n")
	print(x$tile.list)
}
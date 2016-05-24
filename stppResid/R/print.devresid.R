print.devresid <- function(x, ...)
{
	print(x$X)
	cat("Residuals:\n")
	print(x$residuals)
	cat("Spatial grid\n")
	x$grid
}
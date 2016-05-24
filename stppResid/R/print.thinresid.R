print.thinresid <- function(x, ...)
{
	print(x$X)
	cat("Thinning rate:\n")
	print(x$k)
	cat("Thinned residuals:\n")
	print(x$residuals)
	cat("Deleted points:\n")	
	print(x$deleted)
}

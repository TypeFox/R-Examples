print.supresid <- function(x, ...)
{
	print(x$X)
	cat("Superposition rate:\n")
	print(x$k)
	cat("Superposed residuals:\n")
	print(x$residuals)
	cat("Superposed points:\n")
	print(x$super)
}
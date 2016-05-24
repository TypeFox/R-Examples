

print.Daim.vector <- function(x, digits=max(3, getOption("digits") - 3), ...){
	class(x) <- "data.frame"
	print(x, ...)
}
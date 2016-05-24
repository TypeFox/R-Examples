printHG <-
function(x)
{
	cat("Call:\n")
	print(x$call)
	cat("\n")
	
	stat=x[[1]]

	
	op=getOption("scipen")
	options(scipen=7)
	
	cat("HG Test\n")
	print(stat$hyper)
	cat("\n")
	options(scipen=op)
}

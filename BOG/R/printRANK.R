printRANK <-
function(x)
{
	cat("Call:\n")
	print(x$call)
	cat("\n")
	
    stat=x[[1]]
	
	op=getOption("scipen")
	options(scipen=7)
	
	cat("MW Test\n")
	print(stat$rank)
	cat("\n")
	options(scipen=op)			
			

}

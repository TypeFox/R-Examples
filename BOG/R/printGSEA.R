printGSEA <-
function(x)
{
	cat("Call:\n")
	print(x$call)
	cat("\n")
	
	stat=x[[1]]

	
	op=getOption("scipen")
	options(scipen=7)
	
	cat("GSEA Test\n")
	if(is.null(stat$gsea)) cat("Null GSEA")
	else if(!is.null(stat$gsea)){

		print(stat$gsea)
		cat("\n")
	}
	options(scipen=op)
}

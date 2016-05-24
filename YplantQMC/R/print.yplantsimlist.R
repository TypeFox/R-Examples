#'@method print yplantsimlist
#'@S3method print yplantsimlist
print.yplantsimlist <- function(x,...){

	nplants <- length(x)
	cat("Yplant - batch simulation result ('yplantsimlist' object).\n\n")
	cat(paste(c(rep("-",30),"\n"),collapse=""))

	cat("Number of plants :", nplants, "\n\n")
	
	cat("To calculate totals by plant, use the summary() function.\n")
}
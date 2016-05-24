
print.insilico <- function(x,...){
	cat("InSilicoVA fitted object:\n")
	cat(paste(length(x$id), "death processed\n"))
	cat(paste(x$Nsim, "iterations performed, with first", 
			  x$burnin, "iterations discarded\n",
			  trunc((x$Nsim - x$burnin)/x$thin), "iterations saved after thinning\n"))
		if(x$useProbbase){
			cat("Fitted with InterVA4 conditional probability matrix\n")
		}else if(x$keepProbbase.level){
			cat("Fitted with re-estimated InterVA4 conditional probability level table\n")	
		}else{
			cat("Fitted with re-estimating InterVA4 conditional probability matrix\n")	
		}    
}

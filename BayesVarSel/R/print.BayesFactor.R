print.BayesFactor <- function(x, ...){
	cat("---------\n")
	cat("Models:\n")
	print(x$models)
	cat("---------\n")
	cat(paste("Bayes factors (expressed in relation to ",names(x$models)[x$nullmodel],")\n", sep=""))
	print(x$BFi0)
	cat("---------\n")
	cat("Posterior probabilities:\n")
	print(round(x$PostProbi,3))    
  }
    

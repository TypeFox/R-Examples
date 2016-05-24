print.sparseBC <-
function(x,...){	
	cat("Object of class \"sparseBC\"\n")
	cat("Contains the cluster labels for the observations and the features.  Also contains the estimated mean matrix.\n")
	print(data.frame("Number of row clusters" = length(unique(x$Cs)),"Number of column clusters" = length(unique(x$Ds)), "Number of non-zero bicluster means" = length(which(as.vector(abs(x$Mus))>1e-5)),...))
		
	cat("Call:\n\t"); print(x$cl,...)
	invisible(x)
}

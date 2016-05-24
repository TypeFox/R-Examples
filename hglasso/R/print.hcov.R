print.hcov <-
function(x,...){	
	cat("Object of class \"hcov\"\n")
	cat("Contains the estimated matrices V, Z, and Sigma:\n")
	print(data.frame("lambda1" = signif(x$lambda1,2),"lambda2" = signif(x$lambda2,2),"lambda3" = signif(x$lambda3,2), "Number of Edges" = (sum(abs(x$Sigma)!=0)-nrow(x$Sigma))/2, "Indices for hub nodes" = toString(x$hubind)),...)
		
	cat("Call:\n\t");print(x$cl,...)
	invisible(x)
}

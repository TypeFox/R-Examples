print.hglasso <-
function(x,...){	
	cat("Object of class \"hglasso\"\n")
	cat("Contains the estimated matrices V, Z, and Theta:\n")
	print(data.frame("lambda1" = signif(x$lambda1,2),"lambda2" = signif(x$lambda2,2),"lambda3" = signif(x$lambda3,2), "Number of Edges" = (sum(abs(x$Theta)!=0)-nrow(x$Theta))/2, "Indices for hub nodes" = toString(x$hubind)),...)
		
	cat("Call:\n\t");print(x$cl,...)
	invisible(x)
}

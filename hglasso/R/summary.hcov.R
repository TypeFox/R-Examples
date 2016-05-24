summary.hcov <- function(object, ...){
	x <- object
	cat("Summary for the object \"hcov\"\n")
	cat("Information of the optimization problem: \n")
	print(data.frame("n" = x$n, "p" = x$p, "lambda1" = signif(x$lambda1,2),"lambda2" = signif(x$lambda2,2),"lambda3" = signif(x$lambda3,2)),...)
	
	cat("Sigma\n:")
	print(data.frame("Number of Edges" = (sum(abs(x$Sigma)!=0)-x$p)/2, "Indices for hub nodes" = toString(x$hubind)),...)
	cat("V\n:")
	tempV<-(x$V!=0)
	diag(tempV)<-0
	print(data.frame("Indices for hub nodes" = x$hubind, "Number of Edges within each hub" = apply(tempV,2,sum)[x$hubind]),...)
		
	cat("Z\n:")
	print(data.frame("Number of Edges" = (sum(abs(x$Z)!=0)-x$p)/2),...)

	invisible(tempV)	
	invisible(x)
}

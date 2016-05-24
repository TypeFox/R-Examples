print.TCA <-
function(x, ...){
	if(x$method == "pearson")
		cat("Method: Pearson \n")
	if(x$method == "kendall")
		cat("Method: Kendall \n")
	if(x$method == "spearman")
		cat("Method: Spearman \n")
	if(x$method == "npn")
		cat("Method: Nonparanormal (NPN) \n")
	if(x$method == "ns")
		cat("Method: Normal Score (NS) \n")		

	if(x$algorithm == "tp")
		cat("Algorithm: Truncated Power Algorithm \n")
	if(x$algorithm == "spca")
		cat("Algorithm: SPCA algorithm \n")
	if(x$algorithm == "pmd")
		cat("Algorithm: Penalized Matrix Decomposition \n")

	if(x$cov.input) cat("Input: The Covariance Matrix\n")
	if(!x$cov.input) cat("Input: The Data Matrix\n")	
	cat(x$K, "sparse PCs","\n")
	cat("Proportion of variance explained by first k components:", x$pev, "\n")
	cat("Number of nonzero loadings:", apply((x$loadings!=0),2,sum),"\n")
	cat("Sparse Loadings \n")
	print(x$loadings)  
}

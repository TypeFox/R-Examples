print.TCE <-
function(x,...){
   cat("Transelliptical Correlation Estimation via Thresholding \n")
   cat("Method:", x$method, "\n")
	 if(x$cov.input) cat("Input: The Covariance Matrix\n")
	 if(!x$cov.input) cat("Input: The Data Matrix\n")
	
	 cat("Path length:",length(x$lambda),"\n")
	 cat("Sparsity level:",min(x$sparsity),"----->",max(x$sparsity),"\n")    
}

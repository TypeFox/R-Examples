print.summary.blca <-
function(x, ...){
	cat("__________________\n")
	cat("\nBayes-LCA\n")
	cat("Diagnostic Summary\n")
	cat("__________________\n")
	cat("\nHyper-Parameters: \n")
	cat("\n Item Probabilities:\n")
	cat("\n alpha: \n")
	print(x$sum2$ItemProb$alpha)
	cat("\n beta: \n")
	print(x$sum2$ItemProb$beta)
	cat("\n Class Probabilities:\n")
	cat("\n delta: \n")
	print(x$sum2$ClassProb)
	cat("__________________\n")
	cat("\nMethod:",x$method, " \n")
	for(ind in 1:length(x$sum1)){
		cat("\n", x$printnames[ind], x$sum1[ind], "\n")
		}
}

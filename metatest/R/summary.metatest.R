summary.metatest <-
function(object, digits=4, ...) {
	
	variance <- object$variance
	colnames(variance) <- c("parameter", "LLR", "df", "p(LLR)")
	rownames(variance)=c("variance") 
 	
	cat("\nCall: ") 
	print(object$call)
	
 	cat("\n")
 	cat("Number of iterations untill convergence:", object$iter)
 	
 	cat("\n\n")
 	cat("Estimate of residual variance")
 	cat("\n")
 	print(round(variance,digits))
	
	resmat <- object$coefficients
	resmat <- cbind(resmat,object$se,object$tval, object$pZtest, object$pttest, object$dfttest, 
		object$LLR, object$pLLR, object$bartLLR, object$pBartlett, object$ppermtest)
	
	colnames(resmat) <- c("coef","se","t", "p(z-test)", "p(t-test)", "df(t)", "LLR","p(LLR)","B-LLR", "p(B-LLR)", "p(perm)")
	
	cat("\n\n")
	print(round(resmat,digits))	
	cat("\n")
	
	return(invisible(list(var=variance,coef=resmat)))
	
}

print.metatest <-
function(x, ...) {
 	summary(x, ...)
}


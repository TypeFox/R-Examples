summary.blca.em <-
function(object, ...){
	sum1<- c(object$iter, object$poststore[length(object$poststore)] - object$poststore[length(object$poststore)-1], object$logpost, object$AIC, object$BIC)
	names(sum1)<- c("IterNumber", "ConvergenceDiff", "Log-Posterior", "AIC", "BIC")
	
	object$method<- "EM algorithm"
	object$printnames<- c("Number of iterations:","Log-Posterior Increase at Convergence:", "Log-Posterior:", "AIC:", "BIC:")
	object$sum1<- sum1

	NextMethod("summary")
	}

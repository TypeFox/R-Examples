summary.blca.vb <-
function(object, ...){
	sum1<- c(object$iter, object$lbstore[length(object$lbstore)] - object$lbstore[length(object$lbstore)-1], object$LB)
	names(sum1)<- c("IterNumber", "ConvergenceDiff", "Lower Bound")
	
	object$method<- "Variational Bayes"
	object$printnames<- c("Number of iterations:", "Lower Bound Increase at Convergence:", "Lower Bound:")
	object$sum1<- sum1

	NextMethod("summary")
	}

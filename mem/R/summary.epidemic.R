summary.epidemic <-
function(object, ...){
	cat("Call:\n")
	print(object$call)
 	cat("\nOptimum:\n")
	print(object$optimum.map[1:3])
	cat("\nTiming:\n")
	print(object$optimum.map[4:5])
	cat("\nPre-epidemic values:\n")
	print(object$pre.epi)
	cat("\nPost-epidemic values:\n")
	print(object$post.epi)
}

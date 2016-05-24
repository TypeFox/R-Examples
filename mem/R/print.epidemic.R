print.epidemic <-
function(x, ...){
	cat("Call:\n")
	print(x$call)
 	cat("\nOptimum:\n")
	print(x$optimum.map[1])
	cat("\nTiming:\n")
	print(x$optimum.map[4:5])
}

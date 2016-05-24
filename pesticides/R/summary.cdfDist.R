summary.cdfDist <-
function(object, ...){
	cat("Distance:", round(object$cdfDist, 4), "\n")
	cat("No. values distance computed over:", length(object$x), "\n")
}


summary.clustab <-
function(object, ...) {
	x <- object
	if (!inherits(x, "clustab")) 
        	stop("use only with \"clustab\" objects")
	cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  	cat("\n")
}


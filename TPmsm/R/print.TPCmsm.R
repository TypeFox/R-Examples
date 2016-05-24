print.TPCmsm <- function(x, ...) {
	if ( !inherits(x, "TPCmsm") ) stop("'x' must be of class 'TPCmsm'")
	lst <- TransMatrix(x)
	if (x$method == "IPCW2") cat("Inverse Probability of Censoring Weighted conditional transition probabilities\n")
	else if (x$method == "LIN2") cat("Lin conditional transition probabilities\n")
#	cat("\n")
#	cat( paste("Bandwidth = ", x$h, "\n", sep = "") )
	for ( i in 1:length(x$x) ) {
		cat("\n")
		cat( paste("Estimates of P(", x$s, ", ", x$t, " | ", x$x[i], ")\n", sep = "") )
		print(lst[[i]][[1]])
		if ( !is.null(x$n.boot) ) {
			cat("\n")
			cat( paste("Bootstrap confidence bands with", x$n.boot, "samples\n", sep=" ") )
			cat("\n")
			cat( paste( (1-x$conf.level)*50, "%\n", sep="" ) )
			print(lst[[i]][[2]])
			cat("\n")
			cat( paste( (1+x$conf.level)*50, "%\n", sep="" ) )
			print(lst[[i]][[3]])
		}
	}
}

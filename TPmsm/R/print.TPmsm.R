print.TPmsm <- function(x, ...) {
	if ( !inherits(x, "TPmsm") ) stop("'x' must be of class 'TPmsm'")
	lst <- TransMatrix(x)
	if (x$method == "AJ") cat("Aalen-Johansen transition probabilities\n\n")
	else if (x$method == "PAJ") cat("Presmoothed Aalen-Johansen transition probabilities\n\n")
	else if ( x$method %in% c("KMW1", "KMW2") ) cat("Kaplan-Meier Weighted transition probabilities\n\n")
	else if ( x$method %in% c("KMPW1", "KMPW2") ) cat("Presmoothed Kaplan-Meier Weighted transition probabilities\n\n")
	else if (x$method == "IPCW1") cat("Inverse Probability of Censoring Weighted transition probabilities\n\n")
	else if (x$method == "LIN1") cat("Lin transition probabilities\n\n")
	else if (x$method == "LS") {
		cat("Location-Scale transition probabilities\n\n")
#		cat( paste("Bandwidth = ", x$h, "\n\n", sep = "") )
	}
	cat( paste("Estimates of P(", x$s, ", ", x$t, ")\n", sep = "") )
	print(lst[[1]])
	if ( !is.null(x$n.boot) ) {
		cat("\n")
		cat( paste("Bootstrap confidence bands with", x$n.boot, "samples\n", sep=" ") )
		cat("\n")
		cat( paste( (1-x$conf.level)*50, "%\n", sep="" ) )
		print(lst[[2]])
		cat("\n")
		cat( paste( (1+x$conf.level)*50, "%\n", sep="" ) )
		print(lst[[3]])
	}
}

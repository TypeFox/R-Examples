## ===========================
## Methods for stsstatd objects
## ===========================

print.diss.rep <- function(x, ...) {
	criterion <- attr(x,"criterion")
	n <- attr(x,"n")
	quality <- attr(x,"Quality")

	cat(" [>] criterion:",criterion,"\n")
	cat(" [>]", n,"objects in the original data set\n")
	cat(" [>]", length(x),"representative(s)\n")
	cat(" [>] overall quality:", round(quality*100,2),"\n")
	cat(" [>] representative(s) index(es):", x[1:length(x)],"\n")
}

summary.diss.rep <- function(object, ...) {
	criterion <- attr(object,"criterion")
	n <- attr(object,"n")

	cat(" [>] criterion:",criterion,"\n")
	cat(" [>]", n,"objects in the original data set\n")
	cat(" [>]", nrow(object),"representative object(s)\n")
	cat(" [>] statistics for the representative set:\n\n")
	print(attr(object,"Statistics"), digits=3, ...)
	cat("\n    na: number of assigned objects\n")
	cat("    nb: number of objects in the neighborhood\n")
	cat("    SD: sum of the na distances to the representative\n")
	cat("    MD: mean of the na distances to the representative\n")
	cat("    DC: sum of the na distances to the center of the complete set\n")
	cat("    V: discrepancy of the subset\n")
	cat("    Q: quality of the representative\n")
}

## ===========================
## Methods for stsstatd objects
## ===========================

print.stslist.rep <- function(x, ...) {
	criterion <- attr(x,"criterion")
	nbseq <- attr(x,"nbseq")
	quality <- attr(x,"Quality")

	cat("\n [>] criterion:",criterion,"\n")
	cat(" [>]", nbseq,"sequence(s) in the original data set\n")
	cat(" [>]", nrow(x),"representative sequence(s)\n")
	cat(" [>] overall quality:", round(quality*100,2),"\n\n")
	NextMethod(x,...)
}

summary.stslist.rep <- function(object, ...) {
	criterion <- attr(object,"criterion")
	nbseq <- attr(object,"nbseq")
	quality <- attr(object,"Quality")

	cat("\n [>] criterion:",criterion,"\n")
	cat(" [>]", nbseq,"sequence(s) in the original data set\n")
	cat(" [>]", nrow(object),"representative sequences\n")
	cat(" [>] overall quality:", quality,"\n")
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

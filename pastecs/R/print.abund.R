"print.abund" <-
function(x, ...) {
	cat("\nSorting of descriptors according to abundance for:", x$data, "\n\n")
	cat("Coefficient f:", x$f, "\n")
	if (!is.null(x$n)) {						# How many variables do we keep?
		cat("Extraction of: ", x$n, " variable(s) from a total of ", length(x$vr), "\n", sep="")
	} else {
		cat(length(x$vr), " variables sorted\n", sep="")
	}
	cat("\n")
	invisible(x)
}

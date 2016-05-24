"print.escouf" <-
function(x, ...) {
	cat("\nEscoufier's method of equivalent vectors for:", x$data, "\n\n")
	cat("Calculation level:", x$calc.level, "\n")
	cat(x$vars[2], "variable(s) calculated on a total of", x$vars[1], "\n")
	if (!is.null(x$level)) {
		# How many variables do we keep at this level?
		nvars <- length(x$RV[x$RV<x$level])
		cat("Extraction level : ", x$level, " = ", nvars, " variable(s)\n", sep="")
	}
	cat("\n")
	invisible(x)
}

`print.corres` <-
function(x, ...) {

	if (!is(x, "corres")) stop("argument should be a correspondence object")

	cat(paste("\nCall:\ncorres(", x@data$inputName, ")\n\n", sep=""))
    cat("\nEigenvalues (trivial first eigenvalue removed):\n\n")
    cat("    ", x@data$eigenvals, "\n\n")
    cat("Eigenvalue rate, in thousandths:\n\n")
    cat("    ", x@data$eigenrates, "\n\n")

	# print(x@data, quote=F)
	invisible(x)
}


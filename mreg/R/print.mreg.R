"print.mreg" <-
function(x, digits = max(3, getOption("digits") - 3), ...) {
cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(x$coefficients)) {
        cat("Coefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
	if( !is.null(x$root.dispersion)){
	cat("\n\nNuisance Parameters:\n")
	print.default( format(x$nuisance$estimate, digits=digits),print.gap=2,
		quote=FALSE)
	}
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)


}


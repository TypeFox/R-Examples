print.rasch <-
function (x, digits = max(3, getOption("digits") - 3), ...) {
    if (!inherits(x, "rasch"))
        stop("Use only with 'rasch' objects.\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    if (is.matrix(coefs <- x$coef)) {
        coefs <- if (x$IRT.param) IRT.parm(x)$parms else c(coefs[, 1], coefs[1, 2])
        p <- length(coefs) - 1
        if (!x$IRT.param)
            names(coefs) <- c(names(coefs)[1:p], "z")
        cat("Coefficients:\n")
        print(round(coefs, 3), print.gap = 2, quote = FALSE)
    } else
        cat("No coefficients\n")
    cat("\nLog.Lik:", round(x$log.Lik, 3))
    cat("\n\n")
    invisible(x)
}

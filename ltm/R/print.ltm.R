print.ltm <-
function (x, digits = max(3, getOption("digits") - 4), ...) {
    if (!inherits(x, "ltm"))
        stop("Use only with 'ltm' objects.\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    if (is.matrix(coefs <- coef(x))) {
        cat("Coefficients:\n")
        if (x$IRT.param)
            coefs <- IRT.parm(x)$parms
        print(round(coefs, digits), print.gap = 2, quote = FALSE)
    } else
        cat("No coefficients\n")
    cat("\nLog.Lik:", round(x$log.Lik, 3))
    cat("\n\n")
    invisible(x)
}

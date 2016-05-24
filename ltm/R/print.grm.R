print.grm <-
function (x, digits = max(3, getOption("digits") - 4), ...) {
    if (!inherits(x, "grm"))
        stop("Use only with 'grm' objects.\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    coefs <- if (x$IRT.param) IRT.parm(x, digits.abbrv = x$control$digits.abbrv)$parms else x$coef
    coefs <- lapply(coefs, round, digits = digits)
    if (all(sapply(coefs, length) == length(coefs[[1]])))
        coefs <- do.call(rbind, coefs)
    cat("Coefficients:\n")
    print(coefs, print.gap = 2, quote = FALSE)
    cat("\nLog.Lik:", round(x$log.Lik, digits))
    cat("\n\n")
    invisible(x)
}

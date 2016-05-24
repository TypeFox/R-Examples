print.gpcm <-
function (x, digits = max(3, getOption("digits") - 4), ...) {
    if (!inherits(x, "gpcm"))
        stop("Use only with 'gpcm' objects.\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    coefs <- lapply(x$coefficients, round, digits = digits)
    if (all(sapply(coefs, length) == length(coefs[[1]])))
        coefs <- do.call(rbind, coefs)
    cat("Coefficients:\n")
    print(coefs, print.gap = 2, quote = FALSE)
    cat("\nLog.Lik:", round(x$log.Lik, digits))
    cat("\n\n")
    invisible(x)
}

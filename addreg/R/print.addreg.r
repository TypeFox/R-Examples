print.addreg <- function (x, digits = max(3L, getOption("digits") - 3L), ...) 
{
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients")
        if (is.character(co <- x$contrasts)) 
            cat("  [contrasts: ", apply(cbind(names(co), co), 
                1L, paste, collapse = "="), "]")
        cat(":\n")
        print.default(format(x$coefficients, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    if (substr(x$family$family,1,7) == "negbin1") cat("\nScale:", format(x$scale, digits = digits),"\n")
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", 
        x$df.residual, "Residual\n")
    if (nzchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
    cat("Null Deviance:\t   ", format(signif(x$null.deviance, 
        digits)), "\nResidual Deviance:", format(signif(x$deviance, 
        digits)), "\tAIC:", format(signif(x$aic, digits)))
    cat("\n")
    invisible(x)
}
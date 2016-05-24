`print.cusp` <-
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall: ", deparse(x$call), "\n\n")
    if (length(coef(x))) {
        cat("Coefficients")
        if (is.character(co <- x$contrasts)) 
            cat("  [contrasts: ", apply(cbind(names(co), co), 
                1, paste, collapse = "="), "]")
        cat(":\n")
        print.default(format(x$coefficients, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", 
        x$df.residual, "Residual\n")
    if (nchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
    cat("Null Deviance:\t   ", format(signif(x$null.deviance, 
        digits)), paste("\n", colnames(x$fitted), sep = ""), 
        "Deviance:\t", format(signif(x$deviance, digits)), "\tAIC:", 
        format(signif(x$aic, digits)), "\n")
    invisible(x)
}


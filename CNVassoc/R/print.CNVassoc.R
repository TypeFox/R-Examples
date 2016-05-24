print.CNVassoc<-
function (x, digits = 6, ...)
{
    cat("\nCall: ", deparse(x$call), "\n\n")
    cat("Coefficients")
    cat(":\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2,
        quote = FALSE)
    cat("\nNumber of individuals:", length(x$y), "\n")
    cat("Number of estimated parameters:", logLik(x)[2], "\n")
    cat("Deviance:", format(signif(-2 * logLik(x)[1], digits)), "\n")
}
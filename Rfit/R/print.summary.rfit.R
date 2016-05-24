print.summary.rfit <- function (x, digits = max(5, .Options$digits - 2), ...) {
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    est <- x$coef
    printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
    cat("\nMultiple R-squared (Robust):", x$R2, "\n")
    cat("Reduction in Dispersion Test:", round(x$dropstat, digits = digits), 
        "p-value:", round(x$droppval, digits = digits), "\n")
    cat("\n")
}

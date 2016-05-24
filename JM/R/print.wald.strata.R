print.wald.strata <-
function (x, ...) {
    cat("\n\tWald Test for Stratification Factors\n\n")
    pval <- x$Result[, "Pr(> X^2)"]
    pval <- if (pval < 0.0001) "<0.0001" else formatC(round(pval, 4), format = "fg")
    cat("X^2 = ", round(x$Result[, "X^2"], 4), ", df = ", x$Result[, "df"], 
        ", p-value = ", pval, "\n", sep = "")
    cat("alternative hypothesis:", x$alternative, "\n\n")
    invisible(x)
}

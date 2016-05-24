### summaries.R: print and summary functions
### $Id: summaries.R 32 2006-05-17 17:25:23Z bhm $

print.ffmanova <- function(x, digits = max(getOption("digits") - 3, 3), ...) {
    cat("--- 50-50 MANOVA ",
        sum(x$df), " objects -- ",
        ncol(x$pRaw), " responses",
        if (x$stand) " (standardised)",
        ":\n", sep = "")
    tab <- with(x, data.frame(df, exVarSS, c(nPC, NA), c(nBU, NA),
                              c(exVarPC, NA), c(exVarBU, NA), c(pValues, NA)))
    dimnames(tab) <- list(c(x$termNames, "Residuals"),
                          c("Df", "exVarSS", "nPC", "nBU", "exVarPC",
                            "exVarBU", "p-Value"))
    tab <- tab[-1,]                     # Drop the (Intercept) row
    printCoefmat(tab, digits = digits, cs.ind = 2, tst.ind = NULL,
                 zap.ind = c(1,3,4), has.Pvalue = TRUE, na.print = "", ...)
    invisible(x)
}

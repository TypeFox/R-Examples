`print.summary.roc` <- function(x, digits = min(3, getOption("digits") - 4),
                              print.gap = 2, ...) {
    cat("\n")
    writeLines(strwrap("ROC curves of dissimilarities:"))
    cat("\n")
    class(x) <- "data.frame"
    x$`p-value` <- format.pval(x$`p-value`)
    print(x, digits = digits, print.gap = print.gap, ...)
    invisible(x)
}

`print.roc` <- function(x, digits = min(3, getOption("digits") - 4), ...) {
    cat("\n")
    writeLines(strwrap("ROC curve of dissimilarities", prefix = "\t"))
    cat("\n")
    writeLines(strwrap("Discrimination for all groups:"))
    cat("\n")
    with(x$roc$Combined, cat(paste("Optimal Dissimilarity =",
                                  round(optimal, digits), "\n\n")))
    with(x$roc$Combined, cat(paste("AUC = ", round(AUC, digits),
                                  ", p-value: ", format.pval(p.value),
                                  "\n", sep = "")))
    with(x$roc$Combined, cat(paste("No. within:", n.in,
                                  "  No. outside:", n.out, "\n")))
    cat("\n")
    invisible(x)
}

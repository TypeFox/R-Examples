print.descript <-
function (x, digits = max(4, getOption("digits") - 4), ...) {
    if (!inherits(x, "descript"))
        stop("Use only with 'descript' objects.\n")
    cat("\nDescriptive statistics for the", paste("'", x$name, "'", sep = ""), "data-set\n")
    cat("\nSample:\n", x$sample[1], "items and", x$sample[2], "sample units;", if (!is.null(x$missin))
        sum(x$missin[1, ]) else 0, "missing values\n")
    cat("\nProportions for each level of response:\n")
    if (is.list(x$perc))
        print(lapply(x$perc, round, digits = digits))
    else
        print(round(x$perc, digits = digits))
    if (!is.null(x$missin)) {
        cat("\nMissing responses:\n")
        print(round(x$missin, digits = digits))
    }
    cat("\n\nFrequencies of total scores:\n")
    print(x$items)
    if (!is.null(x$bisCorr)) {
        cat("\n\nPoint Biserial correlation with Total Score:\n")
        mat <- round(cbind(x$bisCorr, x$ExBisCorr), digits) 
        colnames(mat) <- c("Included", "Excluded")
        print(mat)
    }
    cat("\n\nCronbach's alpha:\n")
        print(round(x$alpha, digits))
    if (!is.null(x$pw.ass)) {
        cat("\n\nPairwise Associations:\n")
        print(x$pw.ass[seq(1, min(x$n.print, nrow(x$pw.ass))), ])
    }
    cat("\n\n")
    invisible(x)
}

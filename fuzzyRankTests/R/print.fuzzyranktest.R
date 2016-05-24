print.fuzzyranktest <- function(x, digits = 4, ...)
{
    cat("\n")
    writeLines(strwrap(x$method, prefix = "\t"))
    cat("\n")
    cat("data: ", x$data.name, "\n")
    # cat("\n")
    cat("alternative hypothesis: ")
    alt.char <- switch(x$alternative, two.sided = "not equal to",
        less = "less than", greater = "greater than")
    cat("true", names(x$null.value), "is", alt.char, x$null.value, "\n")
    cat("\n")

    if (is.null(x$alpha)) {
        writeLines(strwrap("fuzzy P-value has continuous, piecewise linear CDF with knots and values"))
        cat("\n")
        foobaz <- cbind(x$knots, x$values)
        dimnames(foobaz) <- list(rep("", nrow(foobaz)), c("knots", "values"))
        print.default(foobaz, digits = digits, quote = FALSE, right = TRUE)
    } else {
        writeLines(strwrap(paste("randomized test rejects at level", x$alpha, "with probability", round(x$reject, digits))))
    }
    cat("\n")

    invisible(x)
}

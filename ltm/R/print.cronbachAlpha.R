print.cronbachAlpha <-
function (x, digits = 3, ...) {
    if (!inherits(x, "cronbachAlpha"))
        stop("Use only with 'cronbachAlpha' objects.\n")
    if (x$standardized)
        cat("\nStandardized Cronbach's alpha for the", paste("'", x$name, "'", sep = ""), "data-set\n")
    else
        cat("\nCronbach's alpha for the", paste("'", x$name, "'", sep = ""), "data-set\n")
    cat("\nItems:", x$p)
    cat("\nSample units:", x$n)
    cat("\nalpha:", round(x$alpha, digits))
    if (!is.null(x$ci)) {
        cat("\n\nBootstrap ", 100 * diff(x$probs), "% CI based on ", x$B, " samples\n", sep = "")
        print(round(x$ci, digits))
    }
    cat("\n\n")
    invisible(x)
}

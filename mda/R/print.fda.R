print.fda <-
function (x, ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    cat("\nDimension:", format(x$dimension), "\n")
    cat("\nPercent Between-Group Variance Explained:\n")
    print(round(x$percent, 2))
    error <- x$confusion
    df <- x$fit
    if (!is.null(df)) 
        df <- df$df
    if (!is.null(df)) {
        cat("\nDegrees of Freedom (per dimension):", format(sum(df)), 
            "\n")
    }
    if (!is.null(error)) {
        n <- as.integer(sum(error))
        error <- format(round(attr(error, "error"), 5))
        cat("\nTraining Misclassification Error:", error, "( N =", 
            n, ")\n")
    }
    invisible(x)
}


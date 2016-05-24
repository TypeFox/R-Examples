## User function to compute weighted average tolerances for taxa

## `tolerance` <- function(x, ...)
##    UseMethod("tolerance")

`tolerance.default` <- function(x, env, useN2 = TRUE, ...) {
    x <- data.matrix(x)
    opt <- optima(x, env, ...)
    tol <- sqrt(colSums(x * outer(env, opt, "-")^2) / colSums(x))
    if(useN2) {
        N2 <- sppN2(x)
        tol <- tol / sqrt(1 - (1 / N2))
    }
    names(tol) <- colnames(x)
    class(tol) <- "tolerance"
    attr(tol, "env") <- deparse(substitute(env))
    attr(tol, "N2") <- useN2
    tol
}

`print.tolerance` <- function(x, ...) {
    cat("\n")
    msg <- paste("Weighted Average Tolerances For:", attr(x, "env"))
    writeLines(strwrap(msg, prefix = "\t"),
               sep = "\n\n")
    attr(x, "env") <- NULL
    attr(x, "N2") <- NULL
    print(unclass(x), ...)
}

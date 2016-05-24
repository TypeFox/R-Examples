## simple wrappers to convert optima and tolerance results to
## data frames

`as.data.frame.optima` <- function(x, row.names = NULL,
                                   ...) {
    if (is.matrix(x)) {
        res <- as.data.frame(x)
    } else {
        if(is.null(row.names))
            row.names <- names(x)
        res <- data.frame(Opt = as.numeric(x), row.names = row.names,
                          ...)
    }
    res
}

`as.data.frame.tolerance` <- function(x, row.names = NULL,
                                   ...) {
    if(is.null(row.names))
        row.names <- names(x)
    res <- data.frame(Tol = as.numeric(x), row.names = row.names,
                      ...)
    res
}

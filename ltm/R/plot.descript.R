plot.descript <-
function (x, items = NULL, includeFirstLast = FALSE, xlab, ylab, ...) {
    if (!inherits(x, "descript"))
        stop("Use only with 'descript' objects.\n")
    levs <- apply(data.matrix(x$data)[complete.cases(x$data), ], 2, function(xx) length(unique(xx)))
    if (any(levs > 2))
        stop("the plot method for 'descript' objects currently works for dichotomous responses.\n")
    tot <- as.vector(rowSums(x$data, na.rm = TRUE))
    lis <- split(as.data.frame(x$data), tot)
    if (!includeFirstLast)
        lis <- lis[-c(1, length(lis))]
    out <- sapply(lis, colMeans, na.rm = TRUE)
    if (!is.null(items)) {
        p <- nrow(out)
        if (any(items < 1 | items > p))
            stop("bad value for 'items' argument.\n")
        out <- out[items, ]
    }
    if (missing(xlab))
        xlab <- "Total Score"
    if (missing(ylab))
        ylab <- "Proportion Correct"
    p <- if (!includeFirstLast) 1:ncol(out) else 0:(ncol(out) - 1)
    matplot(cbind(p), t(out), xlab = xlab, ylab = ylab, ...)
    invisible(out)
}

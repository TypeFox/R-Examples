plot.gradientDist <- function(x, orderBy,
                              flipAxes = FALSE,
                              main = NULL,
                              xlab = NULL,
                              ylab = "Distance along gradient",
                              xlim = NULL, ylim = NULL, ...) {
    X <- as.numeric(x)
    if(missing(orderBy)) {
        orderBy <- seq_along(X)
        if(is.null(xlab))
            xlab <- "Sample"
    } else {
        if(is.null(xlab))
            xlab <- deparse(substitute(orderBy))
    }
    xlim <- if(is.null(xlim))
        range(orderBy[is.finite(orderBy)])
    else xlim
    ylim <- if(is.null(ylim))
        range(X[is.finite(X)])
    else ylim
    if(flipAxes)
        plot.default(x = X, y = orderBy, xlab = ylab, ylab = xlab,
                     main = main, ylim = xlim, xlim = ylim, ...)
    else
        plot.default(x = orderBy, y = X, xlab = xlab, ylab = ylab,
                     main = main, ylim = ylim, xlim = xlim, ...)
    invisible(x)
}

lines.gradientDist <- function(x, orderBy, flipAxes = FALSE,
                               type = "l", ...) {
    X <- as.numeric(x)
    if(missing(orderBy)) {
        orderBy <- seq_along(X)
    }
    if(flipAxes)
        lines.default(x = X, y = orderBy, type = type, ...)
    else
        lines.default(x = orderBy, y = X, type = type, ...)
}

points.gradientDist <- function(x, orderBy, flipAxes = FALSE, type = "p",
                               ...) {
    X <- as.numeric(x)
    if(missing(orderBy)) {
        orderBy <- seq_along(X)
    }
    if(flipAxes)
        points.default(x = X, y = orderBy, type = type, ...)
    else
        points.default(x = orderBy, y = X, type = type, ...)
}

##
##  e z p l o t . R
##


ezplot <- function(f, a, b, n = 101, col = "blue", add = FALSE,
                   lty = 1, lwd = 1, marker = 0, pch = 1,
                   grid = TRUE, gridcol = "gray", 
                   fill = FALSE, fillcol = "lightgray",
                   xlab = "x", ylab = "f (x)", main = "Function Plot", ...) {
    fun <- match.fun(f)
    f <- function(x) fun(x)
    stopifnot(is.numeric(a), is.numeric(b),
              length(a) == 1, length(b) == 1, a < b)

    x <- seq(a, b, length.out = n)
    y <- f(x)

    if (!add)
        plot(x, y, type = "n",
        xlab = xlab, ylab = ylab, main = main, ...)

    if (grid)
        grid(col = gridcol)

    if (fill) {
        xx <- c(x, rev(x))
        yy <- c(rep(0, length(x)), rev(y))
        polygon(xx, yy, col = fillcol, border = "darkgray")
    }

    lines(x, y, col = col, lty = lty, lwd = lwd)

    if (marker > 0) {
        m <- min(max(marker, 3), n %/% 3)
        d  <- c(0, sqrt(diff(x)^2 + diff(y)^2))
        cs <- cumsum(d)
        s  <- cs[n]  # sum(d)
        l  <- s / (m-1)

        inds <- numeric(m)
        inds[c(1, m)] <- c(1, n)
        for (k in 2:(m-1))
            inds[k] <- which.min(abs(cs - (k-1)*l))
        points(x[inds], y[inds], col = col, pch = pch)
    } 

    invisible(NULL)
}


ezcontour <- function(f, xlim = c(-pi,pi), ylim = c(-pi,pi), 
                         n = 60, filled = FALSE, col = NULL) {
    fun <- match.fun(f)
    f <- function(x) fun(x)
    stopifnot(is.numeric(xlim), is.numeric(ylim),
              length(xlim) == 2, length(ylim) == 2,
              xlim[1] < xlim[2], ylim[1] < ylim[2])

    xx <- linspace(xlim[1], xlim[2], n)
    yy <- linspace(ylim[1], ylim[2], n)
    F <- matrix(NA, n, n)
    for (i in 1:n) {
        for (j in 1:n) {
            F[i, j] <- f(c(xx[i], yy[j]))
        }
    }
    if (filled) {
        if (is.null(col)) col <- heat.colors(12)
        image(xx, yy, F, col = col)
        contour(xx, yy, F, add = TRUE)
    } else {
        if (is.null(col)) col <- "black"
        contour(xx, yy, F)
        grid()
    }
    invisible(NULL) 
}


ezmesh <- function(f, xlim = c(-pi,pi), ylim = c(-pi,pi), 
                         n = 60, col = "lightgray", ...) {
    fun <- match.fun(f)
    f <- function(x) fun(x)
    stopifnot(is.numeric(xlim), is.numeric(ylim),
              length(xlim) == 2, length(ylim) == 2,
              xlim[1] < xlim[2], ylim[1] < ylim[2])

    xx <- linspace(xlim[1], xlim[2], n)
    yy <- linspace(ylim[1], ylim[2], n)
    F <- matrix(NA, n, n)
    for (i in 1:n) {
        for (j in 1:n) {
            F[i, j] <- f(c(xx[i], yy[j]))
        }
    }
    persp(xx, yy, F, col = col, ...)
    invisible(NULL)
}


ezpolar <- function(fun, interv = c(0, 2*pi)) {
    stopifnot(is.numeric(interv))
    if (length(interv) != 2 || interv[1] >= interv[2])
        stop("Argument 'interv' must have two elements [a, b] with a < b.")

    n <- 91
    x <- seq(interv[1], interv[2], length.out = n)
    y <- fun(x)
    if (length(y) != n) {
        warning("Function 'fun' not vectorized: will do that for you.")
        y <- numeric(n)
        for (i in 1:n) y[i] <- fun(x[i])
    }

    polar(x, y)
    invisible(NULL)
}

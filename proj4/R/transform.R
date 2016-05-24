# Copyright (c) 2007 by Simon Urbanek
#
# part of proj4 R package, license: GPL v2

ptransform <- function(data, src.proj, dst.proj, silent=TRUE) {
    x <- 0
    y <- 0
    z <- 0
    src.proj <- .proj2char(src.proj)
    dst.proj <- .proj2char(dst.proj)
    if (is.list(data)) {
        l <- length(data)
        if (l) {
            x <- data[[1]]
            if (l>1) {
                y <- data[[2]]
                if (l>2) {
                    z <- data[[3]]
                    if (l>3 && !silent) warning("more than three dimensions found, using first three")
                }
            }
        }
    } else {
        d <- dim(data)
        if (is.null(d)) {
            x <- data
        } else {
            if (length(d) != 2 || d[2] > 3)
                stop("input must be 2-dimensional with at most three columns")
            if (d[2] > 0) x <- data[,1]
            if (d[2] > 1) y <- data[,2]
            if (d[2] > 2) z <- data[,3]
        }
    }
    if (!is.numeric(x)) x <- as.numeric(x)
    if (!is.numeric(y)) y <- as.numeric(y)
    if (!is.numeric(z)) z <- as.numeric(z)
    n <- max(length(x),length(y),length(z))
    if (length(x) < n) {
        if (!silent) warning("x is shorter than some of the others, recycling")
        x <- rep(x, length.out=n)
    }
    if (length(y) < n) {
        if (!silent) warning("y is shorter than some of the others, recycling")
        y <- rep(y, length.out=n)
    }
    if (length(z) < n) {
        if (!silent) warning("z is shorter than some of the others, recycling")
        z <- rep(z, length.out=n)
    }
    res <- .C("transform",
              as.character(src.proj), as.character(dst.proj),
              as.integer(n),
              x=as.double(x), y=as.double(y), z=as.double(z),
              NAOK=TRUE, PACKAGE=.package.name)

    if (is.list(data))
        list(x=res$x, y=res$y, z=res$z)
    else
        cbind(res$x, res$y, res$z)
}

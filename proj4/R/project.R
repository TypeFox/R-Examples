# Copyright (c) 2007 by Simon Urbanek
#
# part of proj4 R package, license: GPL v2

project <- function(xy, proj, inverse=FALSE, degrees=TRUE, silent=FALSE, ellps.default="sphere") {
    proj <- .proj2char(proj, ellps.default=ellps.default)
    if (is.list(xy)) {
        if (length(xy)<2) stop("input must be at least 2-dimensional")
        if (length(xy)>2 && !silent) warning("more than two dimensions found, using first two")
        x <- xy[[1]]
        y <- xy[[2]]
    } else {
        d <- dim(xy)
        if (is.null(d) && length(xy)==2) {
            x <- xy[1]
            y <- xy[2]
        } else {
            if (length(d) != 2 || (d[1]!=2 && d[2]!=2))
                stop("input must be 2-dimensional")
            if (d[1]==2) {
                x <- xy[1,]
                y <- xy[2,]
            } else {
                x <- xy[,1]
                y <- xy[,2]
            }
        }
    }
    if (!is.numeric(x)) x <- as.numeric(x)
    if (!is.numeric(y)) y <- as.numeric(y)
    if (length(x) < length(y)) {
        if (!silent) warning("x is shorter than y, recycling")
        x <- rep(x, length.out=length(y))
    }
    if (length(y) < length(x)) {
        if (!silent) warning("y is shorter than x, recycling")
        y <- rep(y, length.out=length(x))
    }
    n <- length(x)
    f <- 0:0
    if (inverse) f <- f + 1:1
    if (degrees) f <- f + 2:2
    res <- .C("project", as.character(proj),                
              as.integer(n), x=as.double(x), y=as.double(y),
              f, NAOK=TRUE, PACKAGE=.package.name)
    if (is.list(xy))
        list(x=res$x, y=res$y)
    else
        cbind(res$x, res$y)
}

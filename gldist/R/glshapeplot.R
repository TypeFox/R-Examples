#  Copyright (C) 2012 Yohan Chalabi
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 or 3 of
#  the License (at your option).
#
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

glshapeplot <- function(x, method, moments = 1:4, ...) {

    #######
    # plot
    plot(c(-1, 1), c(0, 1), type = "n",  ann = FALSE)
    title(main = "GLD Shape Plot", xlab = expression(chi),
          ylab = expression(xi))

    #######
    # add lines of moments
    momentLine <- function(k, lty = 2, col = "grey") {

        small <- 1e-4

        # First part
        chi <- seq( -1 + small, - 2 * sqrt(1 / (4 + k^2)), len = 100)
        b <- chi / (2 * sqrt(1 - chi^2))
        xi <- .5 - sqrt( (1 + 2*b*k + b^2*k^2) / (4 + 8*b*k + k^2 + 4*b^2*k^2) )
        lines(c(-1, chi), c(0, xi), lty = lty, lwd = 2, col = col)

        # second part
        chi <- seq(small - 2 * sqrt(1 / (4 + k^2)), 0, len = 100)
        b <- chi / (2 * sqrt(1 - chi^2))
        xi <- .5 + sqrt( (1 + 2*b*k + b^2*k^2) / (4 + 8*b*k + k^2 + 4*b^2*k^2) )
        lines(chi, xi, lty = lty, lwd = 2, col = col)

        # third part
        chi <- seq(small, 2 * sqrt(1 / (4 + k^2)) - small, len = 100)
        b <- chi / (2 * sqrt(1 - chi^2))
        xi <- .5 + sqrt( (1 - 2*b*k + b^2*k^2) / (4 - 8*b*k + k^2 + 4*b^2*k^2) )
        lines(chi, xi, lty = lty, lwd = 2, col = col)

        # fourth part
        chi <- seq(2 * sqrt(1 / (4 + k^2)), 1 - small)
        b <- chi / (2 * sqrt(1 - chi^2))
        xi <- .5 - sqrt( (1 - 2*b*k + b^2*k^2) / (4 - 8*b*k + k^2 + 4*b^2*k^2) )
        lines(c(chi, 1), c(xi, 0), lty = lty, lwd = 2, col = col)

    }

    if (!is.null(moments))
        sapply(moments, function(m) momentLine(m))

    ans <- NULL

    if (!is.null(x)) {
        if (!is.matrix(x)) x <- as.matrix(x)
        nc <- ncol(x)
        ans <- vector("list", nc)
        for (j in seq.int(nc)) {
            ans[[j]] <- fitgl(x[, j], ...)
            points(ans[[j]]$par[3], ans[[j]]$par[4], pch = 19, cex = .8)
        }
    }

    invisible(ans)
}

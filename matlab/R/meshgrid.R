###
### $Id: meshgrid.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Generate X and Y matrices for three-dimensional plots.
###


##-----------------------------------------------------------------------------
meshgrid <- function(x, y, z, nargout = 2) {

    meshgrid.2d <- function(x, y) {
        if (missing(y)) {
            y <- x
        }

        if (matlab::isempty(x) || matlab::isempty(y)) {
            xx <- matlab::zeros(0, 0)
            yy <- matlab::zeros(0, 0)
        } else {
            xx <- matrix(x, ncol = length(x), byrow = TRUE)
            yy <- matrix(y, nrow = length(y))
            nx <- ncol(xx)
            ny <- nrow(yy)
            xx <- xx[matlab::ones(ny, 1), ]
            yy <- yy[, matlab::ones(1, nx)]
        }

        return(list(x = xx,
                    y = yy))
    }

    meshgrid.3d <- function(x, y, z) {
        if (missing(y) && missing(z)) {
            y <- x
            z <- x
        } else if (missing(z)) {
            stop("not enough input arguments")
        }

        if (matlab::isempty(x) || matlab::isempty(y)) {
            xx <- matlab::zeros(0, 0)
            yy <- matlab::zeros(0, 0)
            zz <- matlab::zeros(0, 0)
        } else {
            nx <- matlab::numel(x)
            ny <- matlab::numel(y)
            nz <- matlab::numel(z)
            xx <- matlab::reshape(as.matrix(x), c(1, nx, 1))
            yy <- matlab::reshape(as.matrix(y), c(ny, 1, 1))
            zz <- matlab::reshape(as.matrix(z), c(1, 1, nz))
            xx <- xx[matlab::ones(ny, 1), , matlab::ones(nz, 1)]
            yy <- yy[, matlab::ones(1, nx), matlab::ones(nz, 1)]
            zz <- zz[matlab::ones(ny, 1), matlab::ones(nx, 1), ]
        }

        return(list(x = xx,
                    y = yy,
                    z = zz))
    }

    if (!is.numeric(nargout)) {
        stop(sprintf("argument %s must be numeric", sQuote("nargout")))
    } else if (!(length(nargout) == 1)) {
        stop(sprintf("argument %s must be of length 1", sQuote("nargout")))
    }

    return(switch(EXPR = as.character(nargout),
                  "2" = meshgrid.2d(x, y),
                  "3" = meshgrid.3d(x, y, z),
                  stop(sprintf("argument %s must be either 2 or 3",
                               sQuote("nargout")))))
}


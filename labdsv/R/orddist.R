orddist <- function (x, dim)
{
    if (class(x) == "pca") {
        z <- x$scores
    } else if (inherits(x, c("pco", "nmds", "metaMDS"))) {
        z <- x$points
    } else if (class(x) == "fso") {
        z <- as.matrix(x$mu)
    }

    if (missing(dim)) dim <- ncol(z)
    if (dim != ncol(z))
        cat(paste("Only comparing first",dim,"dimensions\n"))
    if (dim > ncol(z)) {
        dim <- ncol(z)
        cat(paste("The ordination is only",dim,"dimensionsal."))
    }

    tmp <- dist(z[, 1:dim])
    tmp
}


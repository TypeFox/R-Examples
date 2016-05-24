partition <- function(x, dist, ...)
{
    UseMethod("partition")
}

partition.default <- function (x, dist, ...)
{
    if (class(dist) != 'dist') {
        stop("You must supply an object of class dist as the second argument")
    }
    x <- as.numeric(clustify(x))
    out <- list()
    attr(out,"call") <- match.call()
    out$diss <- dist
    out$clustering <- x
    out$silinfo <- silhouette(x,dist)
    attr(out,'class') <- 'partition'
    return(out)
}

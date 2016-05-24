
### ceeboo 2006

cluster.dist <- function(x, beta) {
    if (!inherits(x, "dist"))
        stop("'x' not of class dist")
    storage.mode(x) <- storage.mode(beta) <- "double"
    obj <- .Call(R_cluster_dist, x, beta)
    names(obj) <- attr(x,"Labels")
    obj
}

###


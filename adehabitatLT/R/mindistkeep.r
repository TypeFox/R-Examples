mindistkeep <- function(x, threshold)
{
    if (!inherits(x,"ltraj"))
        stop("x should be of class 'ltraj'")
    foo <- function(y) {
        ul <- 1
        for (i in 2:nrow(y))
            ul[i] <- ifelse(y$dist[i-1]<threshold,ul[i-1],i)
        if (!is.null(attr(y, "infolocs"))) {
            infol <- attr(y, "infolocs")
        }
        z <- y[ul,c("x","y")]
        da <- y$date
        if (!is.null(attr(y, "infolocs"))) {
            infol <- attr(y, "infolocs")
            lt <- as.ltraj(z, da, id=attr(y,"id"), burst=attr(y,"burst"),
                           typeII=attr(x,"typeII"), infolocs=infol)
        } else {
            lt <- as.ltraj(z, da, id=attr(y,"id"), burst=attr(y,"burst"),
                           typeII=attr(x,"typeII"))
        }
        return(lt)
    }
    res <- do.call("c.ltraj", lapply(x, foo))
    class(res) <- c("ltraj","list")
    return(res)
}

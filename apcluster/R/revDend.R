midDend.local <-
    function (x) if (is.null(mp <- attr(x, "midpoint"))) 0 else mp

memberDend.local <-
    function (x) if (is.null(r <- attr(x, "members"))) 1 else r

isLeaf.local <-
    function (x) (is.logical(L <- attr(x, "leaf"))) && L

midCacheDend.local <- function (x)
{
    stopifnot(inherits(x, "dendrogram"))

    setmid <- function(d)
    {
        if (isLeaf.local(d))
            return(d)

        k <- length(d)

        if (k < 1)
            stop("dendrogram node with non-positive #{branches}")

        r <- d
        midS <- 0

        for (j in 1:k)
        {
            r[[j]] <- unclass(setmid(d[[j]]))
            midS <- midS + midDend.local(r[[j]])
        }

        if (k == 2)
            attr(r, "midpoint") <- (memberDend.local(d[[1]]) + midS) / 2
        else
            attr(r, "midpoint") <- midDend.local(d)

        r
    }

    setmid(x)
}

revDend.local <- function (x)
{
    if (isLeaf.local(x))
        return(x)

    k <- length(x)

    if (k < 1)
        stop("dendrogram non-leaf node with non-positive #{branches}")

    r <- x

    for (j in 1:k)
        r[[j]] <- revDend.local(x[[k + 1 - j]])

    midCacheDend.local(r)
}

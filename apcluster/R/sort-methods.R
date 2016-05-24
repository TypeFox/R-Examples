sort.ExClust <- function(x, decreasing=FALSE,
                         sortBy=c("aggExCluster", "size",
                                  "nameExemplar", "noExemplar"), ...)
{
    sortBy <- match.arg(sortBy)

    if (sortBy == "aggExCluster")
    {
        if (all(dim(x@sim) <= 1))
            stop("cannot sort by agglomerative clustering\n",
                 "because similarity matrix not included in object")
        else
            perm <- aggExCluster(x=x)@order
    }
    else if (sortBy == "size")
        perm <- order(sapply(x@clusters, length))
    else if (sortBy == "nameExemplar")
    {
        if (length(names(x@exemplars)) > 0)
            perm <- order(x@exemplars)
        else
            stop("no names available for exemplars")
    }
    else if (sortBy == "noExemplar")
        perm <- order(x@exemplars)

    if (decreasing)
        perm <- rev(perm)

    x@exemplars <- x@exemplars[perm]
    x@clusters <- x@clusters[perm]

    x
}

#setMethod("sort", signature("ExClust"), sort.ExClust, sealed=TRUE)

setMethod("cutree", signature(tree="AggExResult", k="ANY", h="ANY"),
    function(tree, k, h)
    {
        outObj <- new("ExClust")

        if (!missing(k))
        {
            if (!is.finite(k) || floor(k) != ceiling(k))
                stop("'k' is not an integer number")
            else if (k < 1)
                stop("'k' smaller than 1 does not make sense")
            else if (k > tree@maxNoClusters)
                stop("'k' exceeds maximum number of clusters")
        }
        else if (!missing(h))
        {
            if (!is.finite(h))
                stop("'h' must be numeric")
            else if (h < min(tree@height) ||
                     h > max(tree@height))
                stop("'h' exceeds range of values in 'AggExResult' object",
                     "\nthe range is from ", max(tree@height), " (=> 1 cluster) to ",
                     min(tree@height), " (=> ", tree@maxNoClusters,
                     "clusters)")
            else
                k <- max(which(tree@height >= h))
        }
        else if (missing(k))
            stop("provide either 'k' or 'h'")

        outObj@l         <- tree@l
        outObj@exemplars <- tree@exemplars[[k]]
        outObj@clusters  <- tree@clusters[[k]]
        outObj@idx       <- rep(0, outObj@l)

        for (i in 1:length(outObj@clusters))
        {
            outObj@idx[outObj@clusters[[i]]] <- outObj@exemplars[i]

            if (length(names(outObj@clusters[[i]])) > 0)
                names(outObj@idx)[outObj@clusters[[i]]] <-
                    names(outObj@exemplars)[i]
        }

        outObj@sim  <- tree@sim
        outObj@call <- tree@call

        outObj
    }
)

setMethod("cutree", signature(tree="APResult", k="ANY", h="ANY"),
    function(tree, k, h)
    {
        outObj <- new("ExClust")

        outObj@l         <- tree@l
        outObj@exemplars <- tree@exemplars
        outObj@clusters  <- tree@clusters
        outObj@idx       <- tree@idx
        outObj@sim       <- tree@sim
        outObj@call      <- tree@call

        outObj
    }
)

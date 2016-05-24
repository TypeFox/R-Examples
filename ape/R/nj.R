## nj.R (2009-11-23)

##   Neighbor-Joining Tree Estimation

## Copyright 2004-2009 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

nj <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    if (any(is.na(X)))
        stop("missing values are not allowed in the distance matrix\nConsider using njs()")
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C(C_nj, as.double(X), as.integer(N), integer(2*N - 3),
              integer(2*N - 3), double(2*N - 3), NAOK = TRUE)
    obj <- list(edge = cbind(ans[[3]], ans[[4]]), edge.length = ans[[5]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    reorder(obj)
}

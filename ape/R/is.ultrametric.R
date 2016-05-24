## is.ultrametric.R (2012-02-09)

##   Test if a Tree is Ultrametric

## Copyright 2003-2012 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

is.ultrametric <- function(phy, tol = .Machine$double.eps^0.5)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    if (is.null(phy$edge.length))
        stop("the tree has no branch lengths")

    phy <- reorder(phy)
    n <- length(phy$tip.label)
    e1 <- phy$edge[, 1]
    e2 <- phy$edge[, 2]
    EL <- phy$edge.length

    ## xx: vecteur donnant la distance d'un noeud
    ##     ou d'un tip a partir de la racine
    xx <- numeric(n + phy$Nnode)

    ## the following must start at the root and follow the
    ## edges contiguously; so the tree must be either in cladewise
    ## order (or in pruningwise but the for loop must start from
    ## the bottom of the edge matrix)

    for (i in seq_len(length(e1)))
        xx[e2[i]] <- xx[e1[i]] + EL[i]

    isTRUE(all.equal.numeric(var(xx[1:n]), 0, tolerance = tol))
}

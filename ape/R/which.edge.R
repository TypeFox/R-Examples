## which.edge.R (2014-07-15)

##   Identifies Edges of a Tree

## Copyright 2004-2014 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

getMRCA <- function(phy, tip)
### Find the MRCA of the tips given as `tip'
### (see `root.R' for comments on the code)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    if (is.character(tip)) tip <- which(phy$tip.label %in% tip)
    if (length(tip) < 2) return(NULL)
    Ntip <- length(phy$tip.label)
    seq.nod <- .Call(seq_root2tip, phy$edge, Ntip, phy$Nnode)
    sn <- seq.nod[tip]
    MRCA <- Ntip + 1
    i <- 2
    repeat {
        x <- unique(unlist(lapply(sn, "[", i)))
        if (length(x) != 1) break
        MRCA <- x
        i <- i + 1
    }
    MRCA
}

which.edge <- function(phy, group)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    if (is.character(group))
        group <- which(phy$tip.label %in% group)
    if (length(group) == 1)
        return(match(group, phy$edge[, 2]))

    n <- length(phy$tip.label)
    sn <- .Call(seq_root2tip, phy$edge, n, phy$Nnode)[group]
    i <- 2L
    repeat {
        x <- unique(unlist(lapply(sn, "[", i)))
        if (length(x) != 1) break
        i <- i + 1L
    }
    d <- -(1:(i - 1L))
    x <- unique(unlist(lapply(sn, function(x) x[d])))
    match(x, phy$edge[, 2L])
}

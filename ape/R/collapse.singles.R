## collapse.singles.R (2015-06-22)

##    Collapse "Single" Nodes

## Copyright 2015 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

collapse.singles <- function(tree, root.edge = FALSE)
{
    n <- length(tree$tip.label)
    e1 <- tree$edge[, 1]
    e2 <- tree$edge[, 2]

    tab <- tabulate(e1)
    if (all(tab > 1)) return(tree)

    if (is.null(tree$edge.length)) {
        root.edge <- FALSE
        wbl <- FALSE
    } else {
        wbl <- TRUE
        el <- tree$edge.length
    }

    if (root.edge) ROOTEDGE <- 0

    ## start with the root node:
    ROOT <- n + 1L
    while (tab[ROOT] == 1) {
        i <- which(e1 == ROOT)
        ROOT <- e2[i]
        if (wbl) {
            if (root.edge) ROOTEDGE <- ROOTEDGE + el[i]
            el <- el[-i]
        }
        e1 <- e1[-i]
        e2 <- e2[-i]
    }

    singles <- which(tabulate(e1) == 1)
    while (length(singles)) {
        i <- which(e1 == singles[1])
        j <- which(e2 == e1[i])
        e2[j] <- e2[i]
        if (wbl) {
            el[j] <- el[j] + el[i]
            el <- el[-i]
        }
        e1 <- e1[-i]
        e2 <- e2[-i]
        singles <- which(tabulate(e1) == 1)
    }

    Nnode <- length(e1) - n + 1L

    oldnodes <- unique(e1)
    if (!is.null(tree$node.label))
        tree$node.label <- tree$node.label[oldnodes - n]
    newNb <- integer(max(oldnodes))
    newNb[ROOT] <- n + 1L
    sndcol <- e2 > n
    e2[sndcol] <- newNb[e2[sndcol]] <- n + 2:Nnode
    e1 <- newNb[e1]
    tree$edge <- cbind(e1, e2, deparse.level = 0)
    tree$Nnode <- Nnode
    if (wbl) {
        if (root.edge) tree$root.edge <- ROOTEDGE
        tree$edge.length <- el
    }
    tree
}

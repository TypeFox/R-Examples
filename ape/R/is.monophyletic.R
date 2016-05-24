## is.monophyletic.R (2012-03-23)

##   Test Monophyly

## Copyright 2009-2012 Johan Nylander and Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

is.monophyletic <-
    function(phy, tips, reroot = !is.rooted(phy), plot = FALSE, ...)
{
    if (!inherits(phy, "phylo"))
        stop("object 'phy' is not of class 'phylo'")
    n <- length(phy$tip.label)
    if (length(tips) %in% c(1L, n)) return(TRUE)
    ROOT <- n + 1L
    if (is.numeric(tips)) {
        if (any(tips > n))
            stop("incorrect tip#: should not be greater than the number of tips")
        tips <- sort(as.integer(tips))
    }
    if (is.character(tips))
        tips <- which(phy$tip.label %in% tips)

    if (reroot) {
        outgrp <- phy$tip.label[-tips][1]
        phy <- root(phy, outgroup = outgrp, resolve.root = TRUE)
        rerooted <- TRUE
    } else rerooted <- FALSE

    phy <- reorder(phy)

    seq.nod <- .Call(seq_root2tip, phy$edge, n, phy$Nnode)
    sn <- seq.nod[tips]
    newroot <- ROOT
    i <- 2
    repeat {
        x <- unique(unlist(lapply(sn, "[", i)))
        if (length(x) != 1) break
        newroot <- x
        i <- i + 1
    }
    desc <- which(unlist(lapply(seq.nod, function(x) any(x %in% newroot))))
    if (plot) {
        zoom(phy, tips, subtree = FALSE, ...)
        if (rerooted)
            mtext("Input tree arbitrarily rerooted", side = 1, cex = 0.9)
    }
    ## assuming that both vectors are sorted:
    identical(tips, desc)
}

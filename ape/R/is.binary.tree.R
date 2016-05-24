## is.binary.tree.R (2002-09-12) [modified by EP 2005-05-31, 2005-08-18,
##                                2006-10-04, 2009-05-10]

##    Tests whether a given phylogenetic tree is binary

## Copyright 2002 Korbinian Strimmer

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

is.binary.tree <- function(phy)
{
    if (!inherits(phy, "phylo")) stop('object "phy" is not of class "phylo"')
    ## modified by EP so that it works without edge lengths too (2005-05-31):
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    ## modified by EP so that it works with both rooted and unrooted
    ## trees (2005-08-18):
    if (is.rooted(phy)) {
        if (nb.tip - 1 ==  nb.node) return(TRUE)
        else return(FALSE)
    } else {
        if (nb.tip - 2 ==  nb.node) return(TRUE)
        else return(FALSE)
    }
}

## treeOperators.R (2014-07-09)

##   Trees Operators for Running MCMC

## Copyright 2012-2014 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

getIndexEdge <- function(tip, edge)
    ## 'integer(1)' mustn't be substituted by '0L' except if 'DUP = TRUE':
    .C(get_single_index_integer, as.integer(edge[, 2L]),
       as.integer(tip), integer(1L), NAOK = TRUE)[[3L]]

getIndexEdge2 <- function(node, edge)
    .C(get_two_index_integer, as.integer(edge[, 1L]),
       as.integer(node), integer(2L), NAOK = TRUE)[[3L]]

NeighborhoodRearrangement <- function(phy, n, nodeMax, target, THETA, brtimes)
{
    ## pegas is no more needed:
    ## THETA <- theta.tree(phy, 1)$theta
    bt <- c(rep(0, n), brtimes)
    e <- phy$edge # local copy

### i1, i2, and i3 are edge indices
### target, anc, and sister are node indices

    ## i1 <- which(e2 == target)
    i1 <- getIndexEdge(target, e)
    anc <- e[i1, 1L] # the ancestor of 'target'
    i2 <- getIndexEdge2(target, e) # the 2 edges where 'target' is basal
    i3 <- getIndexEdge2(anc, e) # this includes i1, so:
    i3 <- i3[i3 != i1]
    sister <- e[i3, 2L] # the sister-node of 'target'
    sel <- sample.int(2L, 1L)
    i2.move <- i2[sel]
    i2.stay <- i2[-sel]
    phy$edge[i3, 2L] <- child2move <- e[i2.move, 2L]
    child2stay <- e[i2.stay, 2L]
    phy$edge[i2.move, 2L] <- sister

    ## now adjust branch lengths:
    ## adjust the branch length that was subtending 'sister':
    phy$edge.length[i3] <- bt[anc] - bt[child2move]
    ## random age for 'target' between the ones of 'sister' and 'anc':
    agemin <- max(bt[sister], bt[child2stay])
    pmax <- 1 - exp(-THETA * (bt[anc] - agemin))
    p <- runif(1, 0, pmax)
    newage <- -log(1 - p) / THETA + agemin
    ## alternative with average of ages:
    ## newage <- (bt[anc] + max(bt[sister], bt[child2stay]))/2
    phy$edge.length[i1] <- bt[anc] - newage
    phy$edge.length[i2.move] <- newage - bt[sister]
    ## adjust the branch length below the child that has not been moved:
    phy$edge.length[i2.stay] <- newage - bt[child2stay]

    attr(phy, "order") <- NULL
    phy <- reorder(phy)
    newNb <- integer(nodeMax)
    newNb[n + 1L] <- n + 1L
    sndcol <- phy$edge[, 2L] > n
    phy$edge[sndcol, 2L] <- newNb[phy$edge[sndcol, 2L]] <- (n + 2):nodeMax
    phy$edge[, 1L] <- newNb[phy$edge[, 1L]]
    phy
}

TipInterchange <- function(phy, n)
{
    e <- phy$edge
    repeat {
        k <- sample.int(n, size = 2L)
        i1 <- getIndexEdge(k[1], e)
        i2 <- getIndexEdge(k[2], e)
        ## check that the two tips in 'k' are not sisters
        if (e[i1, 1L] != e[i2, 1L]) break
    }
    e[c(i2, i1), 2L] <- k
    phy$edge <- e
    phy
}

EdgeLengthJittering <- function(phy)
### all edge lengths are added to a random value on U[-MIN, MAX]
### (the ultrametric nature of the tree is kept)
{
    z <- range(phy$edge.length)
    MIN <- z[1]
    MAX <- z[2]
    x <- runif(1, -MIN, MAX) # should be OK even if MIN=0
    phy$edge.length <- phy$edge.length + x
    phy
}

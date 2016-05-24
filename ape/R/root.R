## root.R (2015-11-29)

##   Root of Phylogenetic Trees

## Copyright 2004-2015 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

is.rooted <- function(phy)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    if (!is.null(phy$root.edge)) TRUE
    else
        if (tabulate(phy$edge[, 1])[length(phy$tip.label) + 1] > 2)
            FALSE else TRUE
}

unroot <- function(phy)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    N <- dim(phy$edge)[1]
    if (N < 3)
        stop("cannot unroot a tree with less than three edges.")

    ## delete FIRST the root.edge (in case this is sufficient to
    ## unroot the tree, i.e. there is a multichotomy at the root)
    if (!is.null(phy$root.edge)) phy$root.edge <- NULL
    if (!is.rooted(phy)) return(phy)

    n <- length(phy$tip.label)
    ROOT <- n + 1L

### EDGEROOT[1]: the edge to delete
### EDGEROOT[2]: the target where to stick the edge to delete

### If the tree is in pruningwise (or postorder) order, then
### the last two edges are those connected to the root; the node
### situated in phy$edge[N - 2L, 1L] will be the new root...

    ophy <- attr(phy, "order")
    if (!is.null(ophy) && ophy != "cladewise") {
        NEWROOT <- phy$edge[N - 2L, 1L]
        EDGEROOT <- c(N, N - 1L)
    } else {

### ... otherwise, we remove one of the edges coming from
### the root, and eventually adding the branch length to
### the other one also coming from the root.
### In all cases, the node deleted is the 2nd one (numbered
### nb.tip+2 in 'edge'), so we simply need to renumber the
### nodes by adding 1, except the root (this remains the
### origin of the tree).

        EDGEROOT <- which(phy$edge[, 1L] == ROOT)
        NEWROOT <- ROOT + 1L
    }

    ## make sure EDGEROOT is ordered as described above:
    if (phy$edge[EDGEROOT[1L], 2L] != NEWROOT)
        EDGEROOT <- EDGEROOT[2:1]

    phy$edge <- phy$edge[-EDGEROOT[1L], ]

    s <- phy$edge == NEWROOT # renumber the new root
    phy$edge[s] <- ROOT

    s <- phy$edge > NEWROOT # renumber all nodes greater than the new root
    phy$edge[s] <- phy$edge[s] - 1L

    if (!is.null(phy$edge.length)) {
        phy$edge.length[EDGEROOT[2L]] <-
            phy$edge.length[EDGEROOT[2L]] + phy$edge.length[EDGEROOT[1L]]
        phy$edge.length <- phy$edge.length[-EDGEROOT[1L]]
    }

    phy$Nnode <- phy$Nnode - 1L

    if (!is.null(phy$node.label)) {
        if (NEWROOT == n + 2L)
            phy$node.label <- phy$node.label[-1]
        else {
            lbs <- phy$node.label
            tmp <- lbs[NEWROOT - n]
            lbs <- lbs[-c(1, NEWROOT)]
            phy$node.label <- c(tmp, lbs)
        }
    }
    phy
}

root <- function(phy, outgroup, node = NULL,
                 resolve.root = FALSE, interactive = FALSE)
{
    if (!inherits(phy, "phylo"))
        stop('object not of class "phylo"')
    phy <- reorder(phy)
    n <- length(phy$tip.label)
    ROOT <- n + 1L

    if (interactive) {
        node <- identify(phy)$nodes
        cat("You have set resolve.root =", resolve.root, "\n")
    }

    e1 <- phy$edge[, 1L]
    e2 <- phy$edge[, 2L]
    wbl <- !is.null(phy$edge.length)

    if (!is.null(node)) {
        if (node <= n)
            stop("incorrect node#: should be greater than the number of taxa")
        outgroup <- NULL
        newroot <- node
    } else {
        if (is.numeric(outgroup)) {
            if (any(outgroup > n))
                stop("incorrect taxa#: should not be greater than the number of taxa")
        }
        if (is.character(outgroup)) {
            outgroup <- match(outgroup, phy$tip.label)
            if (anyNA(outgroup))
                stop("specified outgroup not in labels of the tree")
        }
        if (length(outgroup) == n) return(phy)
        outgroup <- sort(outgroup) # used below

        ## First check that the outgroup is monophyletic, unless it has only one tip
        if (length(outgroup) > 1) {
            pp <- prop.part(phy)
            ingroup <- (1:n)[-outgroup]
            newroot <- 0L
            for (i in 2:phy$Nnode) {
                if (identical(pp[[i]], ingroup)) {
                    ## inverted with the next if (... (2013-06-16)
                    newroot <- e1[which(e2 == i + n)]
                    break
                }
                if (identical(pp[[i]], outgroup)) {
                    newroot <- i + n
                    break
                }
            }
            if (!newroot)
                stop("the specified outgroup is not monophyletic")
            MRCA.outgroup <- i + n
        } else newroot <- e1[which(e2 == outgroup)]
    }

    N <- Nedge(phy)
    oldNnode <- phy$Nnode

    if (newroot == ROOT) {
        if (!resolve.root) return(phy) # else (resolve.root == TRUE)
        if (length(outgroup) > 1) outgroup <- MRCA.outgroup
        if (!is.null(node))
            stop("ambiguous resolution of the root node: please specify an explicit outgroup")

        k <- which(e1 == ROOT) # find the basal edges
        if (length(k) > 2) {
            i <- which(e2 == outgroup) # outgroup is always of length 1 here
            j <- k[k != i]
            newnod <- oldNnode + n + 1L
            phy$edge[j, 1] <- newnod

            phy$edge <- rbind(c(ROOT, newnod), phy$edge)
            if (wbl) phy$edge.length <- c(0, phy$edge.length)

            phy$Nnode <- phy$Nnode + 1L
        }
    } else {
        phy$root.edge <- NULL # just in case
        Nclade <- tabulate(e1)[ROOT] # degree of the root node

        ## if only 2 edges connect to the root, we have to fuse them:
        fuseRoot <- Nclade == 2

        INV <- logical(N)
        w <- which(e2 == newroot)
        anc <- e1[w]
        i <- w

        nod <- anc

        if (nod != ROOT) {
            INV[w] <- TRUE
            i <- w - 1L
            repeat {
                if (e2[i] == nod) {
                    if (e1[i] == ROOT) break
                    INV[i] <- TRUE
                    nod <- e1[i]
                }
                i <- i - 1L
            }
        }

        ## we keep the edge leading to the old root if needed:
        if (!fuseRoot) INV[i] <- TRUE

        ## bind the other clades...
        for (j in 1:Nclade) {
            ## do we have to fuse the two basal edges?
            if (fuseRoot) {
                k <- which(e1 == ROOT)
                k <- if (k[2] > w) k[2] else k[1]
                phy$edge[k, 1] <- phy$edge[i, 2]
                if (wbl)
                    phy$edge.length[k] <- phy$edge.length[k] + phy$edge.length[i]
            }
        }

        if (fuseRoot) phy$Nnode <- oldNnode - 1L

        phy$edge[INV, ] <- phy$edge[INV, 2:1]
        if (fuseRoot) {
            phy$edge <- phy$edge[-i, ]
            if (wbl) phy$edge.length <- phy$edge.length[-i]
            N <- N - 1L
        }

        if (resolve.root) {
            newnod <- oldNnode + n + 1L
            if (length(outgroup) == 1L) {
                wh <- which(phy$edge[, 2] == outgroup)
                                        #phy$edge[1] <- newnod
                k <- which(phy$edge[, 1] == newroot) # wh should be among k
                phy$edge[k[k != wh], 1] <- newnod
                o <- c((1:N)[-wh], wh)
                phy$edge <- rbind(c(newroot, newnod), phy$edge[o, ])
                if (wbl) phy$edge.length <- c(0, phy$edge.length[o])
            } else {
                wh <- which(phy$edge[, 1] == newroot)
                phy$edge[wh[-1], 1] <- newnod
                s1 <- 1:(wh[2] - 1)
                s2 <- wh[2]:N
                phy$edge <-
                    rbind(phy$edge[s1, ], c(newroot, newnod), phy$edge[s2, ])
                if (wbl)
                    phy$edge.length <- c(phy$edge.length[s1], 0, phy$edge.length[s2])
            }
            phy$Nnode <- phy$Nnode + 1L
        }
    }

    ## The block below renumbers the nodes so that they conform
    ## to the "phylo" format
    newNb <- integer(n + phy$Nnode)
    newNb[newroot] <- n + 1L
    sndcol <- phy$edge[, 2] > n
    newNb[sort(phy$edge[sndcol, 2])] <- n + 2:phy$Nnode
    phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
    phy$edge[, 1] <- newNb[phy$edge[, 1]]

    if (!is.null(phy$node.label)) {
        newNb <- newNb[-(1:n)]
        if (fuseRoot) {
            newNb <- newNb[-1]
            phy$node.label <- phy$node.label[-1]
        }
        phy$node.label <- phy$node.label[order(newNb)]
        if (resolve.root) {
            phy$node.label[is.na(phy$node.label)] <- phy$node.label[1]
            phy$node.label[1] <- "Root"
        }
    }
    attr(phy, "order") <- NULL
    reorder.phylo(phy)
}

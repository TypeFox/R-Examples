## rtree.R (2013-04-02)

##   Generates Trees

## Copyright 2004-2013 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

rtree <- function(n, rooted = TRUE, tip.label = NULL, br = runif, ...)
{
    foo <- function(n, pos) {
        n1 <- sample.int(n - 1, 1, FALSE, NULL)
        n2 <- n - n1
        po2 <- pos + 2*n1 - 1
        edge[c(pos, po2), 1] <<- nod
        nod <<- nod + 1L
        if (n1 > 2) {
            edge[pos, 2] <<- nod
            foo(n1, pos + 1)
        } else if (n1 == 2) {
            edge[c(pos + 1, pos + 2), 1] <<- edge[pos, 2] <<- nod
            nod <<- nod + 1L
        }
        if (n2 > 2) {
            edge[po2, 2] <<- nod
            foo(n2, po2 + 1)
        } else if (n2 == 2) {
            edge[c(po2 + 1, po2 + 2), 1] <<- edge[po2, 2] <<- nod
            nod <<- nod + 1L
        }
    }

    if (n < 2) stop("a tree must have at least 2 tips.")
    nbr <- 2 * n - 3 + rooted
    edge <- matrix(NA, nbr, 2)

    n <- as.integer(n)
    if (n == 2) {
        if (rooted) edge[] <- c(3L, 3L, 1L, 2L)
        else stop("an unrooted tree must have at least 3 tips.")
    } else if (n == 3) {
        edge[] <-
          if (rooted) c(4L, 5L, 5L, 4L, 5L, 1:3)
          else c(4L, 4L, 4L, 1:3)
    } else if (n == 4 && !rooted) {
        edge[] <- c(5L, 6L, 6L, 5L, 5L, 6L, 1:4)
    } else {
        nod <- n + 1L
        if (rooted) { # n > 3
            foo(n, 1)
            ## The following is slightly more efficient than affecting the
            ## tip numbers in foo(): the gain is 0.006 s for n = 1000.
            i <- which(is.na(edge[, 2]))
            edge[i, 2] <- 1:n
        } else { # n > 4
            n1 <- sample.int(n - 2L, 1L)
            if (n1 == n - 2L) {
                n2 <- n3 <- 1L
            } else {
                n2 <- sample.int(n - n1 - 1L, 1L)
                n3 <- n - n1 - n2
            }
            po2 <- 2L * n1
            po3 <- 2L * (n1 + n2) - 1L
            edge[c(1, po2, po3), 1L] <- nod
            nod <- nod + 1L
            if (n1 > 2L) {
                edge[1L, 2L] <- nod
                foo(n1, 2L)
            } else if (n1 == 2L) {
                edge[2:3, 1L] <- edge[1L, 2L] <- nod
                nod <- nod + 1L
            }
            if (n2 > 2L) {
                edge[po2, 2L] <- nod
                foo(n2, po2 + 1L)
            } else if (n2 == 2L) {
                edge[c(po2 + 1L, po2 + 2), 1L] <- edge[po2, 2L] <- nod
                nod <- nod + 1L
            }
            if (n3 > 2) {
                edge[po3, 2L] <- nod
                foo(n3, po3 + 1L)
            } else if (n3 == 2L) {
                edge[c(po3 + 1L, po3 + 2), 1L] <- edge[po3, 2L] <- nod
                ## nod <- nod + 1L
            }
            i <- which(is.na(edge[, 2L]))
            edge[i, 2] <- 1:n
        }
    }
    phy <- list(edge = edge)
    phy$tip.label <-
      if (is.null(tip.label)) paste("t", sample(n), sep = "")
      else sample(tip.label)
    if (!is.null(br)) {
        phy$edge.length <-
            if (is.function(br)) br(nbr, ...) else rep(br, length.out = nbr)
    }
    phy$Nnode <- n - 2L + rooted
    class(phy) <- "phylo"
    attr(phy, "order") <- "cladewise"
    phy
}

rcoal <- function(n, tip.label = NULL, br = "coalescent", ...)
{
    n <- as.integer(n)
    nbr <- 2*n - 2
    edge <- matrix(NA, nbr, 2)
    ## coalescence times by default:
    x <- if (is.character(br)) 2*rexp(n - 1)/(as.double(n:2) * as.double((n - 1):1))
    else if (is.numeric(br)) rep(br, length.out = n - 1) else br(n - 1, ...)
    if (n == 2) {
        edge[] <- c(3L, 3L, 1:2)
        edge.length <- rep(x, 2)
    } else if (n == 3) {
        edge[] <- c(4L, 5L, 5L, 4L, 5L, 1:3)
        edge.length <- c(x[c(2, 1, 1)], sum(x))
    } else {
        edge.length <- numeric(nbr)
        h <- numeric(2*n - 1)
        node.height <- cumsum(x)
        pool <- 1:n
        nextnode <- 2L*n - 1L
        for (i in 1:(n - 1)) {
            y <- sample(pool, size = 2)
            ind <- (i - 1)*2 + 1:2
            edge[ind, 2] <- y
            edge[ind, 1] <- nextnode
            edge.length[ind] <- node.height[i] - h[y]
            h[nextnode] <- node.height[i]
            pool <- c(pool[! pool %in% y], nextnode)
            nextnode <- nextnode - 1L
        }
    }
    phy <- list(edge = edge, edge.length = edge.length)
    if (is.null(tip.label))
        tip.label <- paste("t", 1:n, sep = "")
    phy$tip.label <- sample(tip.label)
    phy$Nnode <- n - 1L
    class(phy) <- "phylo"
    phy <- reorder(phy)
    ## to avoid crossings when converting with as.hclust:
    phy$edge[phy$edge[, 2] <= n, 2] <- 1:n
    phy
}

rmtree <- function(N, n, rooted = TRUE, tip.label = NULL, br = runif, ...)
{
    a <- replicate(N, rtree(n, rooted = rooted, tip.label =  tip.label,
                            br = br, ...), simplify = FALSE)
    class(a) <- "multiPhylo"
    a
}

stree <- function(n, type = "star", tip.label = NULL)
{
    type <- match.arg(type, c("star", "balanced", "left", "right"))
    n <- as.integer(n)
    if (type == "star") {
        N <- n
        m <- 1L
    } else {
        m <- n - 1L
        N <- n + m - 1L
    }
    edge <- matrix(0L, N, 2)

    switch(type, "star" = {
        edge[, 1] <- n + 1L
        edge[, 2] <- 1:n
    }, "balanced" = {
        if (log2(n) %% 1)
            stop("'n' is not a power of 2: cannot make a balanced tree")
        foo <- function(node, size) {
            if (size == 2) {
                edge[c(i, i + 1L), 1L] <<- node
                edge[c(i, i + 1L), 2L] <<- c(nexttip, nexttip + 1L)
                nexttip <<- nexttip + 2L
                i <<- i + 2L
            } else {
                for (k in 1:2) { # do the 2 subclades
                    edge[i, ] <<- c(node, nextnode)
                    nextnode <<- nextnode + 1L
                    i <<- i + 1L
                    foo(nextnode - 1L, size/2)
                }
            }
        }
        i <- 1L
        nexttip <- 1L
        nextnode <- n + 2L
        foo(n + 1L, n)
    }, "left" = {
        edge[c(seq.int(from = 1, to = N - 1, by = 2), N), 2L] <- 1:n
        nodes <- (n + 1L):(n + m)
        edge[seq.int(from = 2, to = N - 1, by = 2), 2L] <- nodes[-1]
        edge[, 1L] <- rep(nodes, each = 2)
    }, "right" = {
        nodes <- (n + 1L):(n + m)
        edge[, 1L] <- c(nodes, rev(nodes))
        edge[, 2L] <- c(nodes[-1], 1:n)
    })

    if (is.null(tip.label))
        tip.label <- paste("t", 1:n, sep = "")
    phy <- list(edge = edge, tip.label = tip.label, Nnode = m)
    class(phy) <- "phylo"
    attr(phy, "order") <- "cladewise"
    phy
}

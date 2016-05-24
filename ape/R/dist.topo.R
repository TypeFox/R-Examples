## dist.topo.R (2015-04-24)

##      Topological Distances, Tree Bipartitions,
##   Consensus Trees, and Bootstrapping Phylogenies

## Copyright 2005-2015 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.getTreesFromDotdotdot <- function(...)
{
    obj <- list(...)
    if (length(obj) == 1 && !inherits(obj[[1]], "phylo")) obj <- obj[[1]]
    obj
}

dist.topo <- function(x, y = NULL, method = "PH85")
{
    if (!is.null(y)) {
        res <- .dist.topo(x, y, method = method)
    } else {
        n <- length(x) # number of trees
        res <- numeric(n * (n - 1) /2)
        nms <- names(x)
        if (is.null(nms)) nms <- paste0("tree", 1:n)
        k <- 0L
        for (i in 1:(n - 1)) {
            for (j in (i + 1):n) {
                k <- k + 1L
                res[k] <- .dist.topo(x[[i]], x[[j]], method = method)
            }
        }
        attr(res, "Size") <- n
        attr(res, "Labels") <- nms
        attr(res, "Diag") <- attr(res, "Upper") <- FALSE
        attr(res, "method") <- method
        class(res) <- "dist"
    }
    res
}

.dist.topo <- function(x, y, method = "PH85")
{
    if (method == "score" && (is.null(x$edge.length) || is.null(y$edge.length)))
        stop("trees must have branch lengths for branch score distance.")
    nx <- length(x$tip.label)
    x <- unroot(x)
    y <- unroot(y)
    bp1 <- .Call(bipartition, x$edge, nx, x$Nnode)
    bp1 <- lapply(bp1, function(xx) sort(x$tip.label[xx]))
    ny <- length(y$tip.label) # fix by Otto Cordero
    ## fix by Tim Wallstrom:
    bp2.tmp <- .Call(bipartition, y$edge, ny, y$Nnode)
    bp2 <- lapply(bp2.tmp, function(xx) sort(y$tip.label[xx]))
    bp2.comp <- lapply(bp2.tmp, function(xx) setdiff(1:ny, xx))
    bp2.comp <- lapply(bp2.comp, function(xx) sort(y$tip.label[xx]))
    ## End
    q1 <- length(bp1)
    q2 <- length(bp2)
    if (method == "PH85") {
        p <- 0
        for (i in 1:q1) {
            for (j in 1:q2) {
                if (identical(bp1[[i]], bp2[[j]]) | identical(bp1[[i]], bp2.comp[[j]])) {
                    p <- p + 1
                    break
                }
            }
        }
        dT <- q1 + q2 - 2 * p # same than:
        ##dT <- if (q1 == q2) 2*(q1 - p) else 2*(min(q1, q2) - p) + abs(q1 - q2)
    }
    if (method == "score") {
        dT <- 0
        found1 <- FALSE
        found2 <- logical(q2)
        found2[1] <- TRUE
        for (i in 2:q1) {
            for (j in 2:q2) {
                if (identical(bp1[[i]], bp2[[j]]) | identical(bp1[[i]], bp2.comp[[j]])) {
                    dT <- dT + (x$edge.length[which(x$edge[, 2] == nx + i)] -
                                y$edge.length[which(y$edge[, 2] == ny + j)])^2
                    found1 <- found2[j] <- TRUE
                    break
                }
            }
            if (found1) found1 <- FALSE
            else dT <- dT + (x$edge.length[which(x$edge[, 2] == nx + i)])^2
        }
        if (!all(found2))
            dT <- dT + sum((y$edge.length[y$edge[, 2] %in% (ny + which(!found2))])^2)
        dT <- sqrt(dT)
    }
    dT
}

.compressTipLabel <- function(x)
{
    ## 'x' is a list of objects of class "phylo" possibly with no class
    if (!is.null(attr(x, "TipLabel"))) return(x)
    ref <- x[[1]]$tip.label
    n <- length(ref)
    if (length(unique(ref)) != n)
        stop("some tip labels are duplicated in tree no. 1")

    ## serious improvement by Joseph W. Brown!
    relabel <- function (y) {
        label <- y$tip.label
        if (!identical(label, ref)) {
            if (length(label) != length(ref))
                stop("one tree has a different number of tips")
            ilab <- match(label, ref)
            if (any(is.na(ilab)))
                stop("one tree has different tip labels")
            ie <- match(1:n, y$edge[, 2])
            y$edge[ie, 2] <- ilab
        }
        y$tip.label <- NULL
        y
    }
    x <- unclass(x) # another killer improvement by Tucson's hackathon (1/2/2013)
    x <- lapply(x, relabel)
    attr(x, "TipLabel") <- ref
    class(x) <- "multiPhylo"
    x
}

prop.part <- function(..., check.labels = TRUE)
{
    obj <- .getTreesFromDotdotdot(...)
    ntree <- length(obj)
    if (ntree == 1) check.labels <- FALSE
    if (check.labels) obj <- .compressTipLabel(obj) # fix by Klaus Schliep (2011-02-21)
    class(obj) <- NULL # fix by Klaus Schliep (2014-03-06)
    for (i in 1:ntree) storage.mode(obj[[i]]$Nnode) <- "integer"
    class(obj) <- "multiPhylo"
    obj <- .uncompressTipLabel(obj) # fix a bug (2010-11-18)
    clades <- .Call(prop_part, obj, ntree, TRUE)
    attr(clades, "number") <- attr(clades, "number")[1:length(clades)]
    attr(clades, "labels") <- obj[[1]]$tip.label
    class(clades) <- "prop.part"
    clades
}

print.prop.part <- function(x, ...)
{
    if (is.null(attr(x, "labels"))) {
        for (i in 1:length(x)) {
            cat("==>", attr(x, "number")[i], "time(s):")
            print(x[[i]], quote = FALSE)
        }
    } else {
        for (i in 1:length(attr(x, "labels")))
          cat(i, ": ", attr(x, "labels")[i], "\n", sep = "")
        cat("\n")
        for (i in 1:length(x)) {
            cat("==>", attr(x, "number")[i], "time(s):")
            print(x[[i]], quote = FALSE)
        }
    }
}

summary.prop.part <- function(object, ...) attr(object, "number")

plot.prop.part <- function(x, barcol = "blue", leftmar = 4, col = "red", ...)
{
    if (is.null(attr(x, "labels")))
      stop("cannot plot this partition object; see ?prop.part for details.")
    L <- length(x)
    n <- length(attr(x, "labels"))
    layout(matrix(1:2, 2, 1), heights = c(1, 3))
    par(mar = c(0.1, leftmar, 0.1, 0.1))
    one2L <- seq_len(L)
    plot(one2L, attr(x, "number"), type = "h", col = barcol, xlim = c(0.5, L),
         xaxs = "i", xlab = "", ylab = "Frequency", xaxt = "n", bty = "n", ...)
    M <- matrix(0L, L, n)
    for (i in one2L) M[i, x[[i]]] <- 1L
    image.default(one2L, 1:n, M, col = c("white", col), xlab = "", ylab = "", yaxt = "n")
    mtext(attr(x, "labels"), side = 2, at = 1:n, las = 1)
}

prop.clades <- function(phy, ..., part = NULL, rooted = FALSE)
{
    if (is.null(part)) {
        obj <- .getTreesFromDotdotdot(...)
        part <- prop.part(obj, check.labels = TRUE)
    }

    ## until ape 3.0-7 it was assumed implicitly that the labels in phy
    ## are in the same order than in 'part' (bug report by Rupert Collins)
    if (!identical(phy$tip.label, attr(part, "labels"))) {
        i <- match(phy$tip.label, attr(part, "labels"))
        j <- match(seq_len(Ntip(phy)), phy$edge[, 2])
        phy$edge[j, 2] <- i
        phy$tip.label <- attr(part, "labels")
    }
    bp <- prop.part(phy)
    if (!rooted) {
        bp <- postprocess.prop.part(bp)
        part <- postprocess.prop.part(part) # fix by Klaus Schliep
        ## actually the above line in not needed if called from boot.phylo()
    }

    n <- numeric(phy$Nnode)
    for (i in seq_along(bp)) {
        for (j in seq_along(part)) {
            ## we rely on the fact the values returned by prop.part are
            ## sorted and without attributes, so identical can be used:
            if (identical(bp[[i]], part[[j]])) {
                n[i] <- attr(part, "number")[j]
                done <-  TRUE
                break
            }
        }
    }
    n
}

boot.phylo <- function(phy, x, FUN, B = 100, block = 1,
                       trees = FALSE, quiet = FALSE, rooted = FALSE)
{
    if (is.null(dim(x)) || length(dim(x)) != 2)
        stop("the data 'x' must have two dimensions (e.g., a matrix or a data frame)")

    boot.tree <- vector("list", B)

    if (!quiet) # suggestion by Alastair Potts
        progbar <- txtProgressBar(style = 3)

    y <- nc <- ncol(x)

    if (block > 1) {
        a <- seq(1, nc - 1, block)
        b <- seq(block, nc, block)
        y <- mapply(":", a, b, SIMPLIFY = FALSE)
    }

    for (i in 1:B) {
        boot.samp <- unlist(sample(y, replace = TRUE))
        boot.tree[[i]] <- FUN(x[, boot.samp])
        if (!quiet) setTxtProgressBar(progbar, i/B)
    }
    if (!quiet) close(progbar)

    ## for (i in 1:B) storage.mode(boot.tree[[i]]$Nnode) <- "integer"
    ## storage.mode(phy$Nnode) <- "integer"

    if (!quiet) cat("Calculating bootstrap values...")

     if (rooted) {
        pp <- prop.part(boot.tree)
        ans <- prop.clades(phy, part = pp, rooted = rooted)
    } else {
        phy <- reorder(phy, "postorder")
        ints <- phy$edge[, 2] > Ntip(phy)
        ans <- countBipartitions(phy, boot.tree)
        ans <- c(B, ans[order(phy$edge[ints, 2])])
    }

    if (trees) {
        class(boot.tree) <- "multiPhylo"
        ans <- list(BP = ans, trees = boot.tree)
    }
    if (!quiet) cat(" done.\n")
    ans
}

### The next function transforms an object of class "prop.part" so
### that the vectors which are identical in terms of split are aggregated.
### For instance if n = 5 tips, 1:2 and 3:5 actually represent the same
### split though they are different clades. The aggregation is done
### arbitrarily. The call to ONEwise() insures that all splits include
### the first tip.
postprocess.prop.part <- function(x)
{
    n <- length(x[[1]])
    N <- length(x)
    w <- attr(x, "number")

    drop <- logical(N)
    V <- numeric(n)
    for (i in 2:(N - 1)) {
        if (drop[i]) next
        A <- x[[i]]
        for (j in (i + 1):N) {
            if (drop[j]) next
            B <- x[[j]]
            if (length(A) + length(B) != n) next
            V[] <- 0L
            V[A] <- 1L
            V[B] <- 1L
            if (all(V == 1L)) {
                drop[j] <- TRUE
                w[i] <- w[i] + w[j]
            }
        }
    }
    if (any(drop)) {
        labels <- attr(x, "labels")
        x <- x[!drop]
        w <- w[!drop]
        attr(x, "number") <- w
        attr(x, "labels") <- labels
        class(x) <- "prop.part"
    }
    ONEwise(x)
}

### This function changes an object of class "prop.part" so that they
### all include the first tip. For instance if n = 5 tips, 3:5 is
### changed to 1:2.
ONEwise <- function(x)
{
    n <- length(x[[1L]])
    v <- 1:n
    for (i in 2:length(x)) {
        y <- x[[i]]
        if (y[1] != 1) x[[i]] <- v[-y]
    }
    x
}

consensus <- function(..., p = 1, check.labels = TRUE)
{
    foo <- function(ic, node) {
        ## ic: index of 'pp'
        ## node: node number in the final tree
        pool <- pp[[ic]]
        if (ic < m) {
            for (j in (ic + 1):m) {
                wh <- match(pp[[j]], pool)
                if (!any(is.na(wh))) {
                    edge[pos, 1] <<- node
                    pool <- pool[-wh]
                    edge[pos, 2] <<- nextnode <<- nextnode + 1L
                    pos <<- pos + 1L
                    foo(j, nextnode)
                }
            }
        }
        size <- length(pool)
        if (size) {
            ind <- pos:(pos + size - 1)
            edge[ind, 1] <<- node
            edge[ind, 2] <<- pool
            pos <<- pos + size
        }
    }

    obj <- .getTreesFromDotdotdot(...)

    if (!is.null(attr(obj, "TipLabel")))
        labels <- attr(obj, "TipLabel")
    else {
        labels <- obj[[1]]$tip.label
        if (check.labels) obj <- .compressTipLabel(obj)
    }
    ntree <- length(obj)
    ## Get all observed partitions and their frequencies:
    pp <- prop.part(obj, check.labels = FALSE)
    ## Drop the partitions whose frequency is less than 'p':
    if (p == 0.5) p <- 0.5000001 # avoid incompatible splits
    pp <- pp[attr(pp, "number") >= p * ntree]
    ## Get the order of the remaining partitions by decreasing size:
    ind <- sort(unlist(lapply(pp, length)), decreasing = TRUE,
                index.return = TRUE)$ix
    pp <- pp[ind]
    n <- length(labels)
    m <- length(pp)
    edge <- matrix(0L, n + m - 1, 2)
    if (m == 1) {
        edge[, 1] <- n + 1L
        edge[, 2] <- 1:n
    } else {
        nextnode <- n + 1L
        pos <- 1L
        foo(1, nextnode)
    }
    structure(list(edge = edge, tip.label = labels,
              Nnode = m), class = "phylo")
}

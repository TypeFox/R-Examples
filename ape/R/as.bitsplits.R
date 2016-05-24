## as.bitsplits.R (2014-12-11)

##   Conversion Among Split Classes

## Copyright 2011-2014 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

as.bitsplits <- function(x) UseMethod("as.bitsplits")

as.bitsplits.prop.part <- function(x)
{
    foo <- function(vect, RAWVECT) {
        res <- RAWVECT
        for (y in vect) {
            i <- ceiling(y/8)
            res[i] <- res[i] | as.raw(2^(8 - ((y - 1) %% 8) - 1))
        }
        res
    }

    N <- length(x) # number of splits
    n <- length(x[[1]]) # number of tips

    nr <- ceiling(n/8)
    mat <- raw(N * nr)
    dim(mat) <- c(nr, N)

    RAWVECT <- raw(nr)

    for (i in 1:N) mat[, i] <- foo(x[[i]], RAWVECT)

    ## add the n trivial splits of size 1... :
    mat.bis <- raw(n * nr)
    dim(mat.bis) <- c(nr, n)
    for (i in 1:n) mat.bis[, i] <- foo(i, RAWVECT)

    ## ... drop the trivial split of size n... :
    mat <- cbind(mat.bis, mat[, -1, drop = FALSE])

    ## ... update the split frequencies... :
    freq <- attr(x, "number")
    freq <- c(rep(freq[1L], n), freq[-1L])

    ## ... and numbers:
    N <- N + n - 1L

    structure(list(matsplit = mat, labels = attr(x, "labels"),
                   freq =  freq), class = "bitsplits")
}

print.bitsplits <- function(x, ...)
{
    cat('Object of class "bitsplits"\n')
    cat('   ', length(x$labels), 'tips\n')
    cat('   ', length(x$freq), 'partitions\n')
}

sort.bitsplits <- function(x, decreasing = FALSE, ...)
{
    o <- order(x$freq, decreasing = decreasing)
    x$matsplit <- x$matsplit[, o]
    x$freq <- x$freq[o]
    x
}

as.prop.part <- function(x) UseMethod("as.prop.part")

as.prop.part.bitsplits <- function(x)
{
    decodeBitsplits <- function(x) {
        f <- function(y) rev(rawToBits(y)) == as.raw(1)
        which(unlist(lapply(x, f)))
    }
    res <- apply(x$matsplit, 2, decodeBitsplits)
    attr(res, "number") <- x$freq
    attr(res, "labels") <- x$labels
    class(res) <- "prop.part"
    res
}

bitsplits <- function(x)
{
    if (inherits(x, "phylo")) {
        x <- reorder(x, "postorder")
        labs <- x$tip.label
        n <- length(labs)
        m <- x$Nnode
        N <- dim(x$edge)[1]
        nr <- ceiling(n/8)
        nc <- N - n # number of internal edges

        o <- .C(bitsplits_phylo, as.integer(n), as.integer(m),
                as.integer(x$edge), as.integer(N), as.integer(nr),
                raw(nr * nc), NAOK = TRUE)[[6]]
        freq <- rep(1L, nc)

    } else {
        if (!inherits(x, "multiPhylo"))
            stop('x is not of class "phylo" or "multiPhylo"')
        x <- .compressTipLabel(x)
        labs <- attr(x, "TipLabel")
        n <- length(labs)
        nr <- ceiling(n/8)
        ans <- .Call(bitsplits_multiPhylo, x, n, nr)
        nc <- ans[[3]]
        o <- ans[[1]][1:(nr * nc)]
        freq <- ans[[2]][1:nc]
    }

    dim(o) <- c(nr, nc)
    structure(list(matsplit = o, labels = labs, freq = freq),
              class = "bitsplits")
}

countBipartitions <- function(phy, X)
{
    n <- Ntip(phy)
    m <- phy$Nnode
    N <- Nedge(phy)

    SPLIT <- bitsplits(phy)
    nr <- nrow(SPLIT$matsplit)
    nc <- ncol(SPLIT$matsplit)
    freq <- rep(0, nc)
    for (tr in X) {
        tr <- ape::reorder.phylo(tr, "postorder")
        e <- tr$edge
        freq <- .C(CountBipartitionsFromTrees, as.integer(n), as.integer(m),
                   as.integer(e), as.integer(N), as.integer(nr), as.integer(nc),
                   as.raw(SPLIT$matsplit), as.double(freq), NAOK = TRUE)[[8]]
    }
    freq
}

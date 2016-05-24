## evonet.R (2012-09-14)

##   Evolutionary Networks

## Copyright 2011-2012 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

evonet <- function(phy, from, to = NULL)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo".')
    if (!is.rooted(phy))
        warning("the tree is unrooted")
    x <- phy

    if (is.null(to)) {
        if (is.data.frame(from))
            from <- as.matrix(from)
        if (!is.matrix(from))
            stop("'from' must be a matrix or a data frame if 'to' is not given")
        if (ncol(from) > 2) {
            warning("'from' has more than two columns: only the first two will be used.")
            ret <- from[, 1:2]
        } else if (ncol(from) < 2) {
            stop("'from' must have at least two columns")
        } else ret <- from
    } else {
        from <- as.vector(from)
        to <- as.vector(to)
        if (length(from) != length(to))
            stop("'from' and 'to' not of the same length after coercing as vectors")
        ret <- cbind(from, to)
    }

    ## check that values are not out of range:
    storage.mode(ret) <- "integer"
    if (any(is.na(ret)))
        stop("some values are NA's after coercing as integers")
    if (any(ret < 0) || any(ret > Ntip(phy) + phy$Nnode))
        stop("some values are out of range")

    x$reticulation <- ret
    class(x) <- c("evonet", "phylo")
    x
}

as.phylo.evonet <- function(x, ...)
{
    x$reticulation <- NULL
    class(x) <- "phylo"
    x
}

plot.evonet <- function(x, col = "blue", lty = 1, lwd = 1, alpha = 0.5,
                        arrows = 0, arrow.type = "classical", ...)
{
    plot.phylo(as.phylo(x), ...)
    edges(x$reticulation[, 1], x$reticulation[, 2],
          col = rgb(t(col2rgb(col)), alpha = 255 * alpha,
          maxColorValue = 255), lty = lty, lwd = lwd,
          arrows = arrows, type = arrow.type)
}

as.networx.evonet <- function(x, weight = NA, ...)
{
    if (any(x$reticulation <= Ntip(x)))
        stop("some tips are involved in reticulations: cannot convert to \"networx\"")
    x <- reorder(x, "postorder")
    ned <- Nedge(x)
    nrt <- nrow(x$reticulation)
    x$edge <- rbind(x$edge, x$reticulation)
    colnames(x$edge) <- c("oldNodes", "newNodes")
    x$reticulation <- NULL
    x$edge.length <- c(x$edge.length, rep(weight, length.out = nrt))
    x$split <- c(1:ned, 1:nrt)
    class(x) <- c("networx", "phylo")
    x
}

as.network.evonet <- function(x, directed = TRUE, ...)
{
    class(x) <- NULL
    x$edge <- rbind(x$edge, x$reticulation)
    as.network.phylo(x, directed = directed, ...)
}

as.igraph.evonet <- function(x, directed = TRUE, use.labels = TRUE, ...)
{
    class(x) <- NULL
    x$edge <- rbind(x$edge, x$reticulation)
    as.igraph.phylo(x, directed = directed, use.labels = use.labels, ...)
}

print.evonet <- function(x, ...)
{
    nr <- nrow(x$reticulation)
    cat("\n    Evolutionary network with", nr, "reticulation")
    if (nr > 1) cat("s")
    cat("\n\n               --- Base tree ---")
    print.phylo(as.phylo(x))
}

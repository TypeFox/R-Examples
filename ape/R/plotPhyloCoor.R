## plotPhyloCoor.R (2013-03-30)

##   Coordinates of a Tree Plot

## Copyright 2008 Damien de Vienne, 2013 Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

plotPhyloCoor <-
    function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL,
              direction = "rightwards", tip.order = NULL, ...)
{
    Ntip <- length(x$tip.label)
    if (Ntip == 1)
        stop("found only one tip in the tree!")
    Nedge <- dim(x$edge)[1]
    if (any(tabulate(x$edge[, 1]) == 1))
        stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles().")
    Nnode <- x$Nnode
    if (is.null(x$edge.length)) use.edge.length <- FALSE
    phyloORclado <- type %in% c("phylogram", "cladogram")
    horizontal <- direction %in% c("rightwards", "leftwards")
    if (phyloORclado) {
        ## changed by KS:
        yy <- numeric(Ntip + Nnode)
        if (!is.null(tip.order)) {
            yy[tip.order] <- 1:length(tip.order)
        } else {
            x <- reorder(x)
            TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
            yy[TIPS] <- 1:Ntip
        }
    }

    xe <- x$edge
    ## first reorder the tree in cladewise order to avoid cophyloplot() hanging:
    ## x <- reorder(reorder(x), order = "pruningwise") ... maybe not needed anymore (EP)
    x <- reorder(x, order = "postorder")
    ereorder <- match(x$edge[, 2], xe[, 2])

    if (phyloORclado) {
        if (is.null(node.pos)) {
            node.pos <- 1
            if (type == "cladogram" && !use.edge.length)
                node.pos <- 2
        }
        if (node.pos == 1)
            yy <- .C("node_height", as.integer(Ntip), as.integer(Nnode),
                as.integer(x$edge[, 1]), as.integer(x$edge[,
                  2]), as.integer(Nedge), as.double(yy),
                PACKAGE = "ape")[[6]]
        else {
            ans <- .C("node_height_clado", as.integer(Ntip),
                as.integer(Nnode), as.integer(x$edge[, 1]), as.integer(x$edge[,
                  2]), as.integer(Nedge), double(Ntip + Nnode),
                as.double(yy), PACKAGE = "ape")
            xx <- ans[[6]] - 1
            yy <- ans[[7]]
        }
        if (!use.edge.length) {
            if (node.pos != 2)
                xx <- .C("node_depth", as.integer(Ntip), as.integer(Nnode),
                         as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
                         as.integer(Nedge), double(Ntip + Nnode), 1L,
                         PACKAGE = "ape")[[6]] - 1
            xx <- max(xx) - xx
        } else {
            xx <- .C("node_depth_edgelength", as.integer(Ntip),
                as.integer(Nnode), as.integer(x$edge[, 1]), as.integer(x$edge[,
                  2]), as.integer(Nedge), as.double(x$edge.length),
                double(Ntip + Nnode), PACKAGE = "ape")[[7]]
        }
    }
    ##if (type == "fan") {
    ##    TIPS <- xe[which(xe[, 2] <= Ntip), 2]
    ##    xx <- seq(0, 2 * pi * (1 - 1/Ntip), 2 * pi/Ntip)
    ##    theta <- double(Ntip)
    ##    theta[TIPS] <- xx
    ##    theta <- c(theta, numeric(Nnode))
    ##    theta <- .C("node_height", as.integer(Ntip), as.integer(Nnode),
    ##        as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
    ##        as.integer(Nedge), theta, DUP = FALSE, PACKAGE = "ape")[[6]]
    ##    if (use.edge.length) {
    ##        r <- .C("node_depth_edgelength", as.integer(Ntip),
    ##            as.integer(Nnode), as.integer(x$edge[, 1]), as.integer(x$edge[,
    ##              2]), as.integer(Nedge), as.double(x$edge.length),
    ##            double(Ntip + Nnode), DUP = FALSE, PACKAGE = "ape")[[7]]
    ##    }
    ##    else {
    ##        r <- .C("node_depth", as.integer(Ntip), as.integer(Nnode),
    ##            as.integer(x$edge[, 1]), as.integer(x$edge[,
    ##              2]), as.integer(Nedge), double(Ntip + Nnode),
    ##            DUP = FALSE, PACKAGE = "ape")[[6]]
    ##        r <- 1/r
    ##    }
    ##    xx <- r * cos(theta)
    ##    yy <- r * sin(theta)
    ##}
    ##if (type == "unrooted") {
    ##    XY <- if (use.edge.length)
    ##        unrooted.xy(Ntip, Nnode, x$edge, x$edge.length)
    ##    else unrooted.xy(Ntip, Nnode, x$edge, rep(1, Nedge))
    ##    xx <- XY$M[, 1] - min(XY$M[, 1])
    ##    yy <- XY$M[, 2] - min(XY$M[, 2])
    ##}
    ##if (type == "radial") {
    ##    X <- .C("node_depth", as.integer(Ntip), as.integer(Nnode),
    ##        as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
    ##        as.integer(Nedge), double(Ntip + Nnode), DUP = FALSE,
    ##        PACKAGE = "ape")[[6]]
    ##    X[X == 1] <- 0
    ##    X <- 1 - X/Ntip
    ##    yy <- c((1:Ntip) * 2 * pi/Ntip, rep(0, Nnode))
    ##    Y <- .C("node_height", as.integer(Ntip), as.integer(Nnode),
    ##        as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
    ##        as.integer(Nedge), as.double(yy), DUP = FALSE, PACKAGE = "ape")[[6]]
    ##    xx <- X * cos(Y)
    ##    yy <- X * sin(Y)
    ##}
    if (phyloORclado && direction != "rightwards") {
        if (direction == "leftwards") {
            xx <- -xx
            xx <- xx - min(xx)
        }
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
            if (direction == "downwards") {
                yy <- -yy
                yy <- yy - min(yy)
            }
        }
    }
    cbind(xx, yy)
}

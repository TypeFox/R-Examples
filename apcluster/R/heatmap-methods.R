heatmap.ExClust <- function(x, y, ...)
{
    if (!all(dim(x@sim) <= 1))
        return(invisible(heatmap(x, x@sim, ...)))
    else
        stop("similarity matrix is missing for heatmap plotting")
}

setMethod("heatmap", signature(x="ExClust", y="missing"), heatmap.ExClust)


heatmap.ExClust.matrix <- function(x, y, ...)
{
    if (all(dim(y) <= 1))
        stop("'y' must be a non-empty similarity matrix")
    else if (x@l != nrow(y))
        stop("size of clustering result does not fit to size of data set")
    else if (nrow(y) != ncol(y) && length(x@sel) == 0)
        stop("similarity matrix must be quadratic")

    if (any(y == -Inf))
    {
        rng <- range(y[which(y > -Inf)])
        fill <- 2 * rng[1] - rng[2]
        y[which(y == -Inf)] <- fill
    }

    aggres <- aggExCluster(y, x)

    heatmap(aggres, y, ...)

    return(invisible(aggres))
}

setMethod("heatmap", signature(x="ExClust", y="matrix"), heatmap.ExClust.matrix)


heatmap.ExClust.sparseMatrix <- function(x, y, ...)
{
    if (all(dim(y) <= 1))
        stop("'y' must be a non-empty similarity matrix")
    else if (x@l != nrow(y))
        stop("size of clustering result does not fit to size of data set")
    else if (nrow(y) != ncol(y))
        stop("similarity matrix must be quadratic")

    aggres <- aggExCluster(y, x, includeSim=TRUE)

    heatmap(aggres, aggres@sim, ...)

    aggres@sim <- matrix(ncol=0, nrow=0)

    return(invisible(aggres))
}

setMethod("heatmap", signature(x="ExClust", y="sparseMatrix"),
          heatmap.ExClust.sparseMatrix)


heatmap.ExClust.Matrix <- function(x, y, ...)
{
    y <- try(as(y, "matrix"))

    if (class(y) == "try-error")
        stop("cannot cast 'y' (class '", class(y), "') to class 'matrix'")

    return(invisible(heatmap(x, y, ...)))
}

setMethod("heatmap", signature(x="ExClust", y="Matrix"), heatmap.ExClust.Matrix)


heatmap.AggExResult <- function(x, y, ...)
{
    if (!all(dim(x@sim) <= 1))
        heatmap(x, x@sim, ...)
    else
        stop("similarity matrix is missing for heatmap plotting")
}

setMethod("heatmap", signature(x="AggExResult", y="missing"),
          heatmap.AggExResult)


heatmap.AggExResult.matrix <- function(x, y, Rowv=TRUE, Colv=TRUE,
                                       sideColors=NULL, col=heat.colors(12),
                                       base=0.05, add.expr, margins=c(5, 5, 2),
                                       cexRow=max(min(35 / nrow(y), 1), 0.1),
                                       cexCol=max(min(35 / ncol(y), 1), 0.1),
                                       main=NULL, dendScale=1, barScale=1,
                                       legend=c("none", "col"),
                                       ...)
{
    if (all(dim(y) <= 1))
        stop("'y' must be a non-empty matrix")
    else if (x@l != nrow(y))
        stop("size of clustering result does not fit to size of data set")
    else if (length(x@sel) == 0 && ncol(y) != nrow(y))
        stop("'y' must be quadratic")
    else if (length(x@sel) > 0 && ncol(y) != length(x@sel))
        stop("no. of columns in 'y' and no. of selected samples in 'x' ",
             "do not match")

    legend <- match.arg(legend)

    rowInd <- unlist(x@clusters[[x@maxNoClusters]][x@order])

    doRdend <- TRUE
    doCdend <- TRUE

    dend <- NULL

    if (is.na(Rowv) || identical(Rowv, FALSE) || x@maxNoClusters < 3)
        doRdend <- FALSE
    else
    {
        dend <- as.dendrogram(x, base=base, useNames=FALSE)
        rowInd <- as.numeric(order.dendrogram(dend))
    }

    colInd <- rowInd

    if (length(x@sel) > 0)
    {
        colInd <- rank(intersect(rowInd, x@sel))
        doCdend <- FALSE
    }
    else if (is.na(Colv) || identical(Colv, FALSE) || x@maxNoClusters < 3)
        doCdend <- FALSE
    else if (!is.na(Colv) && !doRdend)
    {
        dend <- as.dendrogram(x, base=base, useNames=FALSE)
        rowInd <- as.numeric(order.dendrogram(dend))
        colInd <- rowInd
    }

    if ((doRdend || doCdend) && (!is.numeric(dendScale) ||
                                 length(dendScale) != 1 || dendScale <= 0 ||
                                 dendScale > 2))
        stop("'dendScale' must be a single positive value not larger than 2")

    if (is.null(sideColors))
    {
        if (length(x) != nrow(y))
        {
            lx <- length(x)
            lx2 <- lx + if (lx %% 2) 1 else 0
            ind <- as.vector(t(matrix(1:lx2, lx2 / 2)))[1:lx]
            sideColors <- rainbow(lx)[ind]
        }
    }
    else if (any(is.na(sideColors)))
        sideColors <- NULL
    else if (!is.character(sideColors))
        stop("'sideColors' must be vector of colors, NA or NULL")
    else
    {
        if (length(sideColors) < 2)
            stop("use at least two different colors in 'sideColors' argument")

        if (length(sideColors) < length(x))
            sideColors <- rep(sideColors, length.out=length(x))
        else
            sideColors <- sideColors[1:length(x)]
    }

    if (length(rownames(y)) == 0)
    {
        labRow <- as.character(rowInd)

        if (length(x@sel) > 0)
            labCol <- as.character(intersect(rowInd, x@sel))
    }
    else
    {
        labRow <- rownames(y)[rowInd]

        if (length(x@sel) > 0)
            labCol <- rownames(y)[intersect(rowInd, x@sel)]
    }

    if (length(colnames(y)) == 0)
        labCol <- as.character(colInd)
    else
        labCol <- colnames(y)[colInd]

    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) dendScale else 0.05, 4)
    lhei <- c((if (doCdend) dendScale else 0.05) +
              if (!is.null(main)) 0.2 else 0, 4)

    if (length(sideColors) > 0)
    {
        if (!is.numeric(barScale) || length(barScale) != 1 || barScale <= 0 ||
            barScale > 4)
            stop("'barScale' must be a single positive value not larger than 4")

        invIndex <- rep(0, nrow(y))
        for (i in 1:x@maxNoClusters)
            invIndex[x@clusters[[x@maxNoClusters]][[i]]] <- i

        srtIndex <- unique(invIndex[rowInd])

        rowColors <- rep(sideColors,
                         sapply(x@clusters[[x@maxNoClusters]][srtIndex],
                                length))

        if (length(x@sel) > 0)
            colColors <-
                rep(sideColors, sapply(x@clusters[[x@maxNoClusters]][srtIndex],
                                function(cl) length(intersect(cl, x@sel))))
        else
            colColors <- rowColors

        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        lhei <- c(lhei[1L], 0.1 * barScale, lhei[2L])
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1),
                      lmat[, 2] + 1)
        lwid <- c(lwid[1L], 0.1 * barScale, lwid[2L])
    }

    lmat[is.na(lmat)] <- 0

    if (legend != "none")
    {
        lmat <- cbind(lmat, c(rep(0, nrow(lmat) - 1), max(lmat) + 1))
        lwid <- c(lwid, 0.25)
    }

    dev.hold()
    on.exit(dev.flush())
    op <- par(no.readonly=TRUE)
    on.exit(par(op), add=TRUE)
    layout(lmat, widths=lwid, heights=lhei, respect=TRUE)

    if (length(sideColors) > 0)
    {
        par(mar=c(margins[1], 0, 0, 0.5))
        image(rbind(1:nrow(y)), col=rev(rowColors), axes=FALSE)
        par(mar=c(0.5, 0, 0, margins[2]))
        image(cbind(1:ncol(y)), col=colColors, axes=FALSE)
    }

    par(mar=c(margins[1], 0, 0, margins[2]))

    image(1:ncol(y), 1:nrow(y), t(y[rev(rowInd), colInd]),
          xlim=(0.5 + c(0, ncol(y))), ylim=(0.5 + c(0, nrow(y))),
          axes=FALSE, xlab="", ylab="", col=col, ...)

    if (cexCol > 0)
        axis(1, 1:ncol(y), labels=labCol, las=2, line=-0.5, tick=0,
             cex.axis=cexCol)

    if (cexRow > 0)
        axis(4, 1:nrow(y), labels=rev(labRow), las=2, line=-0.5, tick=0,
             cex.axis=cexRow)

    if (!missing(add.expr))
        eval.parent(substitute(add.expr))

    par(mar=c(margins[1], 0, 0, 0))

    if (doRdend)
        plot(revDend.local(dend), horiz=TRUE, axes=FALSE, yaxs="i",
             leaflab="none")
    else
        frame()

    par(mar=c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))

    if (doCdend)
        plot(dend, axes=FALSE, xaxs="i", leaflab="none")
    else
        frame()

    if (!is.null(main))
    {
        par(xpd=NA)
        title(main, cex.main=(1.5 * op[["cex.main"]]))
    }

    if (legend != "none")
    {
        par(mar=c(margins[1], 0, 0, margins[3]))

        rng <- range(y)
        colvals <- seq(rng[1], rng[2], length.out=length(col))

        image(y=colvals, z=rbind(colvals),
              col=col, axes=FALSE, xlab="", ylab="")
        axis(4)
    }

    return(invisible(dend))
}

setMethod("heatmap", signature(x="AggExResult", y="matrix"),
          heatmap.AggExResult.matrix)


setMethod("heatmap", signature(x="matrix", y="missing"),
          function(x, y, ...) stats::heatmap(x=x, ...))

setMethod("heatmap", signature(x="missing", y="matrix"),
          function(x, y, ...) stats::heatmap(x=y, ...))

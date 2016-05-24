ht <-
function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
    distfun = dist, hclustfun = hclust, reorderfun = function(d, 
        w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, 
        "Rowv"), scale = c("row", "column", "none"), na.rm = FALSE, 
    margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 1, cexCol = 1, labRow = NULL, 
    labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
    verbose = getOption("verbose"), q, g) 
{
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (is.null(Rowv)) 
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) 
        Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram")) 
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            if (!is.logical(Rowv) || Rowv) 
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr))) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1L:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram")) 
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc) 
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm) 
                x
            else t(x)))
            ddc <- as.dendrogram(hcc)
            if (!is.logical(Colv) || Colv) 
                ddc <- reorderfun(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc))) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1L:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow)) 
        if (is.null(rownames(x))) 
            (1L:nr)[rowInd]
        else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol)) 
        if (is.null(colnames(x))) 
            (1L:nc)[colInd]
        else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/", check.margin = FALSE)
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/", check.margin = FALSE)
    }
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, 4)
    lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
        4)
    if (!missing(ColSideColors)) {
        if (!is.character(ColSideColors) || length(ColSideColors) != 
            nc) 
            stop("'ColSideColors' must be a character vector of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors) != 
            nr) 
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
            1), lmat[, 2] + 1)
        lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei, 
            "; lmat=\n")
        print(lmat)
    }
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.0))
        image(rbind(1L:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0, 0, 0, margins[2]))
        image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none") 
        x <- t(x)
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1L:nr
    #add extra space to right margin of plot within frame
    par(mar=c(5, 4, 0, 4) + 0.1)
    image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "")
    points(q,g, pch=3, col="black", cex=3, lwd=3)
    axis(1, 1L:nc, labels = labCol, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 3)
    axis(2, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 2, line =+ 1.9)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    par(mar = c(margins[1], 0, 0, 0))
    if (doRdend) 
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
    if (doCdend) 
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else if (!is.null(main)) 
        frame()
    if (!is.null(main)) {
        par(xpd = NA)
        title(main, cex.main = 1)
        #.5 * op[["cex.main"]]
    }
    invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
        doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
} # end ht


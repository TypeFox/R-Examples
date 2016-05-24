heatmap_fn <- function(inputdata, 
    # Dendrograms
    Rowv=TRUE, Colv=if(symm) "Rowv" else TRUE, distfun=dist, hclustfun=hclust,
    dendrogram=c("both", "row", "column", "none"), symm=FALSE,
    # Scaling of data matrix
    scale=c("none", "row", "column"), na.rm=TRUE, 
    # Image plot
    revC=identical(Colv, "Rowv"), 
    # Other plotting-type calls (i.e. points())
    add.expr, 
    # Colour mapping
    breaks, symbreaks=min(inputdata < 0, na.rm=TRUE) || scale!="none", 
    # Colours for heatmap (colramp in HeatMap() function)
    col=NULL, 
    # Separation between colour blocks
    colsep, rowsep, sepcolor="white", sepwidth=c(0.05, 0.05), 
    # Labelling of blocks
    cellnote, notecex=1, notecol="cyan", na.color=par("bg"),
    # Plot level information on colour key
    lvtrace=c("column", "row", "both", "none"), 
    tracecol="cyan", 
    hline=median(breaks), vline=median(breaks), linecol=tracecol, 
    # Row and Column labels
    margins=c(5, 5), 
    # Sample (group), Metabolite colour bars
    ColSideColors, RowSideColors, 
    # Character size etc
    cexRow=0.2 + (1 / log10(numrows)), cexCol=0.2 + (1 / log10(numcols)), 
    labRow=NULL, labCol=NULL, 
    # Colour key, density information
    key=TRUE, keysize=1.5, density.info=c("histogram", "density", "none"),
    denscol=tracecol, symkey=min(inputdata < 0, na.rm=TRUE) || symbreaks, 
    densadj=0.25, 
    # Plot labels
    main=NULL, xlab=NULL, ylab=NULL, 
    # Plot layout
    lmat=NULL, lhei=NULL, lwid=NULL, ...)
{
    # TODO: if no row dendrogram, plot colour scale as strip on LHS.
    
    # Define scaling function
    # @param x:    Data to scale
    # @type x:    data.frame
    scale01 <- function(x, low=min(x), high=max(x)) 
    {
        x <- (x - low) / (high - low)
        invisible(x) ## Just in case it returns
    }
    
    # Define invalid function to replace import from gtools
    gtinvalid <- function (val) 
    {
        if (missing(val) || is.null(val) || length(val) == 0) 
            return(TRUE)
        if (is.list(val)) 
            return(all(sapply(val, gtinvalid)))
        else if (is.vector(val)) 
            return(all(is.na(val)))
        else return(FALSE)
    }
    
    #
    #    Argument matching
    #
    scale <- if (symm && missing(scale)) {
        "none"
    } else {
        match.arg(scale)
    }
    dendrogram <- match.arg(dendrogram)
    lvtrace <- match.arg(lvtrace)
    density.info <- match.arg(density.info)
    tracecol <- if (is.null(tracecol)) {
        "cyan"
    }
    
    #
    #    Sanity checks
    #
    # This should be sufficient - I can't think of what else you may need
    # to test in order to ensure a proper colour palette.
    if (is.null(col)) {
        col=redgreen(75)
    }
    if (is.null(Rowv) || is.na(Rowv)) {
        Rowv <- FALSE
    }
    if (is.null(Colv) || is.na(Colv)) {
        Colv <- FALSE
    } else if (Colv == "Rowv" && !isTRUE(Rowv)) {
        Colv <- FALSE
    }
    if (length(dim(inputdata)) != 2 || !is.numeric(inputdata)) {
        stop("`inputdata' must be a numeric matrix")
    }
    
    if (!missing(breaks) && (scale != "none")) {
        warning(
            paste("Using scale=\"row\" or scale=\"column\" when breaks are",
                "specified can\nproduce unpredictable results.", 
                "Please consider using only one or the other."
            )
        )
    }
    # Check type of inputdata
    numrows <- nrow(inputdata)
    numcols <- ncol(inputdata)
    if (numrows <= 1 || numcols <= 1) {
        stop("`inputdata' must have at least 2 rows and 2 columns")
    }
    if (!is.numeric(margins) || length(margins) != 2) {
        stop("`margins' must be a numeric vector of length 2")
    }
    if (missing(cellnote)) {
        cellnote <- matrix("", ncol=ncol(inputdata), nrow=nrow(inputdata))
    }
    
    # Check for sensible dendrogram plotting arguments
    draw_row_dendro <- dendrogram %in% c("both", "row")
    draw_col_dendro <- dendrogram %in% c("both", "column")
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && draw_row_dendro) {
            # If Colv is a logical, and true
            if (is.logical(Colv) && (Colv)) {
                dendrogram <- "column"
            } else {
                dedrogram <- "none"
            }
            # If row dendrogram requested, while Rowv is false, warn user
            warning(
                paste("Discrepancy: Rowv is FALSE, while dendrogram is `",
                    dendrogram, "'.\nOmitting row dendogram.", sep=""
                )
            )
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && draw_col_dendro) {
            if (is.logical(Rowv) && (Rowv)) {
                dendrogram <- "row"
            } else {
                dendrogram <- "none"
            }
            warning(
                paste("Discrepancy: Colv is FALSE, while dendrogram is `",
                    dendrogram, "'.\nOmitting column dendogram.", sep=""
                )
            )
        }
    }
    
    #
    #    Dendrograms
    #
    # If Rowv is a dendrogram, order the rows by that
    if (inherits(Rowv, "dendrogram")) {
        row_dendro <- Rowv
        rowInd <- order.dendrogram(row_dendro)
    # or if (list of) integers, by that list
    } else if (is.integer(Rowv)) {
        row_hclust <- hclustfun(distfun(inputdata))
        row_dendro <- as.dendrogram(row_hclust)
        row_dendro <- reorder(row_dendro, Rowv)
        rowInd <- order.dendrogram(row_dendro)
        if (numrows != length(rowInd)) {
            stop("Row dendrogram ordering gave index of wrong length")
        }
    # otherwise create dendrogram and reorder
    } else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(inputdata, na.rm=na.rm)
        row_hclust <- hclustfun(distfun(inputdata))
        row_dendro <- reorder(as.dendrogram(row_hclust), Rowv)
        rowInd <- order.dendrogram(row_dendro)
        if (numrows != length(rowInd)) {
            stop("Row dendrogram ordering gave index of wrong length")
        }
    # or just reverse the order of the rows         ## why?
    } else {
        rowInd <- numrows:1
    }
    # Similarly for column dendrogram
    if (inherits(Colv, "dendrogram")) {
        col_dendro <- Colv
        colInd <- order.dendrogram(col_dendro)
    } else if (identical(Colv, "Rowv")) {
        if (numrows != numcols) {
            stop("Colv = \"Rowv\" but nrow(inputdata) != ncol(inputdata)")
        }
        if (exists("row_dendro")) {
            col_dendro <- row_dendro
            colInd <- order.dendrogram(col_dendro)
        } else {
            colInd <- rowInd
        }
    } else if (is.integer(Colv)) {
        if (symm) {
            col_hclust <- hclustfun(distfun(inputdata))
        } else {
            col_hclust <- hclustfun(distfun(t(inputdata)))
        }
        col_dendro <- reorder(as.dendrogram(col_hclust), Colv)
        colInd <- order.dendrogram(col_dendro)
        if (numcols != length(colInd)) {
            stop("Column dendrogram ordering gave index of wrong length")
        }
    } else if (isTRUE(Colv)) {
        Colv <- colMeans(inputdata, na.rm=na.rm)
        if (symm) {
            col_hclust <- hclustfun(distfun(inputdata))
        } else {
            col_hclust <- hclustfun(distfun(t(inputdata)))
        }
        col_dendro <- as.dendrogram(col_hclust)
        col_dendro <- reorder(col_dendro, Colv)
        colInd <- order.dendrogram(col_dendro)
        if (numcols != length(colInd)) {
            stop("Column dendrogram ordering gave index of wrong length")
        }
    } else {
        colInd <- 1:numcols
    }
    
    # Create a data structure to contain all the useful info
    retval <- list()
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()             ## only used for return
    retval$unscaled <- inputdata[rowInd, colInd] # ordered unscaled data
    cellnote <- cellnote[rowInd, colInd]
    # If no row labels specified, create them
    if (is.null(labRow)) {
        labRow <- if (is.null(rownames(inputdata))) {
            (1:numrows)[rowInd]
        } else {
            #rownames(inputdata)[rowInd]
            rownames(retval$unscaled)
        }
    } else {
        labRow <- labRow[rowInd]
    }
    # Similarly for columns
    if (is.null(labCol)) {
        labCol <- if (is.null(colnames(inputdata))) {
            (1:numcols)[colInd]
        } else {
            #colnames(inputdata)[colInd]
            colnames(retval$unscaled)
        }
    } else {
        labCol <- labCol[colInd]
    }
    
    # Data scaling
    if (scale == "row") {
        retval$rowMeans <- rowMeans(retval$unscaled, na.rm=na.rm)
        rowmean_swept_id <- sweep(retval$unscaled, 1, retval$rowMeans)
        retval$rowSDs <- apply(rowmean_swept_id, 1, sd, na.rm=na.rm)
        retval$data_mat <- sweep(rowmean_swept_id, 1, retval$rowSDs, "/")
    } else if (scale == "column") {
        retval$colMeans <- colMeans(retval$unscaled, na.rm=na.rm)
        colmean_swept_id <- sweep(retval$unscaled, 2, retval$colMeans)
        retval$colSDs <- apply(colmean_swept_id, 2, sd, na.rm=na.rm)
        retval$data_mat <- sweep(colmean_swept_id, 2, retval$colSDs, "/")
    } else {
        retval$data_mat <- retval$unscaled
    }
    
    # Determine the number of subdivisions of the data matrix colours
    # (i.e. how many colours to use, and at what values to assign them)
    ## If no 'breaks' in function definition, assign it a value
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
            breaks <- length(col) + 1
    }
    # If 'breaks' is a value (i.e. assigned as such in function definition,
    # or assigned in previous if statement)
    if (length(breaks) == 1) { 
        # If the data is not to be treated symmetrically, define
        # 'breaks' to be between min and max of data, otherwise take the 
        # largest absolute value and make symmetrical about zero.
        if (!symbreaks) {
            breaks <- seq(
                min(retval$data_mat, na.rm=na.rm), 
                max(retval$data_mat, na.rm=na.rm),
                length=breaks
            )
        } else {
            extreme <- max(abs(retval$data_mat), na.rm=TRUE)
            breaks <- seq(-extreme, extreme, length=breaks)
        }
    }
    
    retval$data_mat[retval$data_mat < min(breaks)] <- min(breaks)
    retval$data_mat[retval$data_mat > max(breaks)] <- max(breaks)
    if (exists("row_dendro")) {
        retval$rowDendrogram <- row_dendro
    }
    if (exists("col_dendro")) {
        retval$colDendrogram <- col_dendro
    }
    retval$breaks <- breaks
    retval$col <- col
    
    #
    #    Plot layout
    #
    # Plan layout of plotting device based on what is requested
    if (missing(lhei) || is.null(lhei)) {
        lhei <- c(keysize, 4)
    }
    # TODO: if no row dendro, lwid[1] == 0, lwid[2] == 0.2 to add in scalebar
    #       (see below)
    if (missing(lwid) || is.null(lwid)) {
        lwid <- c(keysize, 4)
    }
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            csc_check <- length(ColSideColors) == numcols
            if (!is.character(ColSideColors) || !csc_check) {
                stop(
                    paste("'ColSideColors' must be a character",
                        "vector of length ncol(inputdata)"
                    )
                )
            }
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            rsc_check <- length(RowSideColors) == numrows
            if (!is.character(RowSideColors) || !rsc_check) {
                stop(
                    paste("'RowSideColors' must be a character",
                        "vector of length nrow(inputdata)"
                    )
                )
            }
            lmat <- cbind(lmat[, 1] + 1, 
                c(rep(NA, nrow(lmat) - 1), 1), 
                lmat[, 2] + 1
            )
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        # Replace any remaining NA with 0 (don't plot there)
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) {
        stop("'lhei' must have length = nrow(lmat) = ", nrow(lmat))
    }
    if (length(lwid) != ncol(lmat)) {
        stop("'lwid' must have length = ncol(lmat) = ", ncol(lmat))
    }
    #op <- par(no.readonly=TRUE) ### Causes errors with par()$pin for some reason
    op <- par(
        cex=par()$cex, cex.axis=par()$cex.axis, cex.main=par()$cex.main, 
        col=par()$col, las=par()$las, lty=par()$lty, lwd=par()$lwd, 
        mar=par()$mar, usr=par()$usr, xaxs=par()$xaxs, yaxs=par()$yaxs, 
        xaxt=par()$xaxt, yaxt=par()$yaxt
    )
    on.exit(par(op))
    
    # layout() - divide device up into subplots
    # layout.show(n) - show the divisions of the subplots
    # plotmap<-layout(lmat, widths=lwid, heights=lhei, respect=FALSE)
    # layout.show(plotmap)
    #
    # TODO: If no row dendro requested:
    #     lwid[,1] <- 0.5 # make a narrow band at left of plot layout
    #     lmat[1,1] <- 0 # don't plot in usual colour key area
    #     #image4zone to be filled with rotated colour key (or lmat[,1])
    # TODO: Option csc/rsc for dendro-side
    #                                  
    # lmat: plotting matrix; lwid: widths of subplots; lhei: heights of same
    # > lmat
    #      [,1] [,2] [,3]
    # [1,]    6    0    5  ## does not plot at 0 points.
    # [2,]    0    0    2  ## centre pos for 3 if csc/rsc?
    # [3,]    4    1    3
    #                      ## order of plots given by order of numbers
    
    layout(lmat, widths=lwid, heights=lhei, respect=FALSE)
    # Plot 1 - Row Side colours
    if (!missing(RowSideColors)) {
        # if (dendrosidecols) {
        # } else {
        par(mar=c(margins[1], 0, 0, 0.5))
        image(rbind(1:numrows), col=RowSideColors[rowInd], axes=FALSE)
        #}
    }
    
    # Plot 2 - Col Side colours
    if (!missing(ColSideColors)) {
        # if dendrosidecols) {
        # } else {
        par(mar=c(0.5, 0, 0, margins[2]))
        image(cbind(1:numcols), col=ColSideColors[colInd], axes=FALSE)
        #}
    }
    
    # Plot 3 - Heatmap
    # TODO: space (label) margins by length of label?
    par(mar=c(margins[1], 0, 0, margins[2]))
    retval$t_data_mat <- t(retval$data_mat)
    cellnote <- t(cellnote)
    if (revC) {
        y_index <- numrows:1
        if (exists("row_dendro")) {
            row_dendro <- rev(row_dendro)
        }
        retval$t_data_mat <- retval$t_data_mat[, y_index]
        cellnote <- cellnote[, y_index]
    } else {
        y_index <- 1:numrows
    }
    # x-pos, y-pos, value to plot, ...
    image(1:numcols, 1:numrows, retval$t_data_mat, 
        xlim=0.5 + c(0, numcols), 
        ylim=0.5 + c(0, numrows), 
        axes=FALSE, xlab="", ylab="", col=col, breaks=breaks, ...
    )
    
    if (!gtinvalid(na.color) & any(is.na(retval$t_data_mat))) {
        mmat <- ifelse(is.na(retval$t_data_mat), 1, NA)
        image(1:numcols, 1:numrows, mmat, axes=FALSE, xlab="", ylab="",
            col=na.color, add=TRUE
        )
    }
    
    # Plot x-axis (sample names)
    # Label axis with (input data) row names (i.e. sample names)
    # May have to change line if csc/rsc on dendroside
    axis(1, 1:numcols, labels=labCol, las=2, line=(-0.5), tick=0, 
        cex.axis=cexCol
    )
    
    # x-axis label
    if (!is.null(xlab)) {
        mtext(xlab, side=1, line=margins[1] - 1.25)
    }
    # Plot y-axis (met names)
    # Label axis with (input data) column names (i.e. metabolite names)
    axis(4, y_index, labels=labRow, las=2, line=(-0.5), tick=0, 
        cex.axis=cexRow
    )
    
    # y-axis label
    if (!is.null(ylab)) {
        mtext(ylab, side=4, line=margins[2] - 1.25)
    }
    
    ### Add other things to plot
    if (!missing(add.expr)) {
        # e.g. if add.expr=points(1:5, 1:5, pch=17, col="yellow")
        #      it will draw points in a diagonal line across the map.
        # NB: the centre of each cell is positioned at integer values.
        eval(substitute(add.expr))
    }
    
    # Column separators (list of columns to draw a line around)
    if (!missing(colsep)) {
        for (csep in colsep) {
            rect(xleft=csep + 0.5, 
                ybottom=rep(0, length(csep)), 
                xright=csep + 0.5 + sepwidth[1], 
                ytop=rep(ncol(retval$t_data_mat) + 1, csep), 
                lty=1, lwd=1, col=sepcolor, border=sepcolor
            )
        }
    }
    # Row separators
    if (!missing(rowsep)) {
        for (rsep in rowsep) {
            rect(xleft=0, 
                ybottom=(ncol(retval$t_data_mat) + 1 - rsep) - 0.5, 
                xright=nrow(retval$t_data_mat) + 1, 
                ytop=(ncol(retval$t_data_mat) + 1 - rsep) - 0.5 - sepwidth[2], 
                lty=1, lwd=1, col=sepcolor, border=sepcolor
            )
        }
    }
    
    #
    #    Level tracing & Cell annotation
    #
    retval$scaled <- scale01(retval$data_mat, min(breaks), max(breaks))
    # Level trace by columns
    if (lvtrace %in% c("both", "column")) {
        retval$vline <- vline                   ## used for return
        vline.vals <- scale01(vline, min(breaks), max(breaks))
        for (ii in colInd) {
            if (!is.null(vline)) {
                abline(v=ii - 0.5 + vline.vals, col=linecol, lty=2)
            }
            # Get x & y trace positions for plotting
            xtpos <- rep(ii, nrow(retval$scaled)) + retval$scaled[, ii] - 0.5
            xtpos <- c(xtpos[1], xtpos)
            ytpos <- 1:length(xtpos) - 0.5
            lines(x=xtpos, y=ytpos, lwd=1, col=tracecol, type="s")
        }
    }
    # Level trace by rows
    if (lvtrace %in% c("both", "row")) {
        retval$hline <- hline                   ## used for return
        hline.vals <- scale01(hline, min(breaks), max(breaks))
        for (ii in rowInd) {
            if (!is.null(hline)) {
                abline(h=ii + hline, col=linecol, lty=2)
            }
            ytpos <- rep(ii, ncol(retval$scaled)) + retval$scaled[ii, ] - 0.5
            ytpos <- rev(c(ytpos[1], ytpos))
            xtpos <- length(ytpos):1 - 0.5
            lines(x=xtpos, y=ytpos, lwd=1, col=tracecol, type="s")
        }
    }
    # Cell labels
    if (!missing(cellnote)) {
        text(x=c(row(cellnote)), y=c(col(cellnote)), labels=c(cellnote),
            col=notecol, cex=notecex
        )
    }
    
    #
    #    Draw dendrograms
    #
    # For metabolites
    par(mar=c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(row_dendro, horiz=TRUE, axes=FALSE, yaxs="i", leaflab="none")
    } else {
        # plot nothing
        plot.new()
    }
    # For samples
    par(mar=c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(col_dendro, axes=FALSE, xaxs="i", leaflab="none")
    } else {
        # plot nothing
        plot.new()
    }
    
    # Plot title
    if (!is.null(main)) {
        title(main, cex.main=1.5 * op[["cex.main"]])
    }
    # Colour key
    # TODO: Determine how this needs to be adjusted for no row dendro
    if (key) {
        par(mar=c(5, 4, 2, 1), cex=0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(retval$t_data_mat, breaks)), na.rm=TRUE)
            min.raw <- -max.raw
            max.mat <- max(abs(retval$t_data_mat), na.rm=TRUE)
            tmpbreaks[1] <- -max.mat
            tmpbreaks[length(tmpbreaks)] <- max.mat
        } else {
            min.raw <- min(retval$t_data_mat, na.rm=TRUE)
            max.raw <- max(retval$t_data_mat, na.rm=TRUE)
        }
        keyvals <- matrix(seq(min.raw, max.raw, length=length(col)), ncol=1)
        # Draw colour ramp
        image(keyvals, col=col, breaks=tmpbreaks, xaxt="n", yaxt="n")
        # Simplify the axis plotting by scaling to 0 -> 1
        par(usr=c(0, 1, 0, 1))
        # Get x-axis labels
        key_xlab <- pretty(breaks)
        # Determine where they should be plotted on the 0 -> 1 scale
        key_x_at <- scale01(as.numeric(key_xlab), min.raw, max.raw)
        # Add x-axis
        axis(1, at=key_x_at, labels=key_xlab)
        
        # Add density plot
        if (density.info == "density") {
            key_dens <- density(retval$t_data_mat, adjust=densadj, na.rm=TRUE)
            omit <- key_dens$x < min(breaks) | key_dens$x > max(breaks)
            key_dens$x <- key_dens$x[-omit]
            key_dens$y <- key_dens$y[-omit]
            key_dens$x <- scale01(key_dens$x, min.raw, max.raw)
            lines(key_dens$x, 
                (key_dens$y / max(key_dens$y)) * 0.95, 
                col=denscol, lwd=1
            )
            axis(2, 
                at=(pretty(key_dens$y) / max(key_dens$y)) * 0.95, 
                pretty(key_dens$y)
            )
            title("Color Key\nand Density Plot")
            par(cex=0.5)
            mtext(side=2, "Density", line=2)
        } else if (density.info == "histogram") {
            key_his <- hist(retval$t_data_mat, plot=FALSE, breaks=breaks)
            # Scale xvals to 0 -> 1
            key_his_xvals <- scale01(breaks, min.raw, max.raw)
            
            key_his_yvals <- c(key_his$counts, 
                key_his$counts[length(key_his$counts)]
            )
            lines(key_his_xvals, 
                (key_his_yvals / max(key_his_yvals)) * 0.95, 
                lwd=1, type="s", col=denscol
            )
            axis(2, 
                at=(pretty(key_his_yvals) / max(key_his_yvals)) * 0.95, 
                pretty(key_his_yvals)
            )
            title("Color Key\nand Histogram")
            par(cex=0.5)
            mtext(side=2, "Count", line=2)
        }
    } else {
        # If no key requested, plot nothing
        plot.new()
    }
    retval$colorTable <- data.frame(
        low=retval$breaks[-length(retval$breaks)],
        high=retval$breaks[-1], color=retval$col
    )
    invisible(retval)
}

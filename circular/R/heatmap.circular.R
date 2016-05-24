#############################################################
#                                                           #
#   heatmap.circular function                               #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: October, 7, 2007                                  #
#   Version: 0.7                                            #
#                                                           #
#   Copyright (C) 2007 Claudio Agostinelli                  #
#                                                           #
#############################################################
## This is a modified version of the heatmap function in package stats.
## The original heatmap function is made by Andy Liaw; modified RG, MM :
#############################################################

heatmap.circular <- function (x, Rowv=NULL, Colv=if(symm)"Rowv" else NULL,
	  distfun = dist.circular, hclustfun = hclust,
          reorderfun = function(d,w) reorder(d,w),
          add.expr, symm = FALSE, revC = identical(Colv, "Rowv"),
	  na.rm=TRUE,
	  margins = c(5, 5), lwid=c(1,4), lhei=c(1,4),
          ColSideColors, RowSideColors,  NAColors="black",
	  cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc),
	  labRow = NULL, labCol = NULL, main = NULL, xlab = NULL, ylab = NULL,
	  keep.dendro = FALSE, 
          annotate.expr, annotate=rep(NA, 4),
	  verbose = getOption("verbose"), ...) {
    if(length(di <- dim(x)) != 2 || !is.numeric(x))
	stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if(nr <= 1 || nc <= 1)
	stop("'x' must have at least 2 rows and 2 columns")
    if(!is.numeric(margins) || length(margins) != 2)
	stop("'margins' must be a numeric vector of length 2")    
    doRdend <- !identical(Rowv,NA)
    doCdend <- !identical(Colv,NA)
    ## by default order by row/col means
    if(is.null(Rowv)) Rowv <- rowMeans(x, na.rm = na.rm)
    if(is.null(Colv)) Colv <- colMeans(x, na.rm = na.rm)

    ## get the dendrograms and reordering indices

    if(doRdend) {
	if(inherits(Rowv, "dendrogram"))
	    ddr <- Rowv
	else {
	    hcr <- hclustfun(distfun(x))
	    ddr <- as.dendrogram(hcr)
	    if(!is.logical(Rowv) || Rowv)
		ddr <- reorderfun(ddr, Rowv)
	}
	if(nr != length(rowInd <- order.dendrogram(ddr)))
	    stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1L:nr

    if(doCdend) {
	if(inherits(Colv, "dendrogram"))
	    ddc <- Colv
	else if(identical(Colv, "Rowv")) {
	    if(nr != nc)
		stop('Colv = "Rowv" but nrow(x) != ncol(x)')
	    ddc <- ddr
	}
	else {
	    hcc <- hclustfun(distfun(if(symm)x else t(x)))
	    ddc <- as.dendrogram(hcc)
	    if(!is.logical(Colv) || Colv)
		ddc <- reorderfun(ddc, Colv)
	}
	if(nc != length(colInd <- order.dendrogram(ddc)))
	    stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1L:nc

    ## reorder x
    x <- x[rowInd, colInd]

    labRow <-
	if(is.null(labRow))
	    if(is.null(rownames(x))) (1L:nr)[rowInd] else rownames(x)
	else labRow[rowInd]
    labCol <-
	if(is.null(labCol))
	    if(is.null(colnames(x))) (1L:nc)[colInd] else colnames(x)
	else labCol[colInd]

    ## Calculate the plot layout
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if(doRdend) lwid[1] else 0.05, lwid[2])
    lhei <- c((if(doCdend) lhei[1] else 0.05) + if(!is.null(main)) 0.2 else 0, lhei[2])
    if(!missing(ColSideColors)) { ## add middle row to layout
	if(!is.character(ColSideColors) || length(ColSideColors) != nc)
	    stop("'ColSideColors' must be a character vector of length ncol(x)")
	lmat <- rbind(lmat[1,]+1, c(NA,1), lmat[2,]+1)
	lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if(!missing(RowSideColors)) { ## add middle column to layout
	if(!is.character(RowSideColors) || length(RowSideColors) != nr)
	    stop("'RowSideColors' must be a character vector of length nrow(x)")
	lmat <- cbind(lmat[,1]+1, c(rep(NA, nrow(lmat)-1), 1), lmat[,2]+1)
	lwid <- c(lwid[1], 0.2, lwid[2])
    }
    # Annotate setup
    if (!is.na(annotate[1])) {
      lmat <- rbind(lmat, c(rep(NA, ncol(lmat)-1), max(lmat, na.rm=TRUE)+1))
      lhei <- c(lhei, annotate[1])
    }
    if (!is.na(annotate[3])) {
      lmat <- rbind(c(rep(NA, ncol(lmat)-1), max(lmat, na.rm=TRUE)+1), lmat)
      lhei <- c(annotate[3], lhei)
    }
    if (!is.na(annotate[2])) {
      lmat <- cbind(c(rep(NA, nrow(lmat)-1-!is.na(annotate[1])), max(lmat, na.rm=TRUE)+1, rep(NA, !is.na(annotate[1]))), lmat)
      lwid <- c(annotate[2], lwid)
    }
    for (i in 4:length(annotate)) {
      if (!is.na(annotate[i])) {
        lmat <- cbind(lmat, c(rep(NA, nrow(lmat)-1-!is.na(annotate[1])), max(lmat, na.rm=TRUE)+1, rep(NA, !is.na(annotate[1]))))
        lwid <- c(lwid, annotate[i])
      }
    }
    
    lmat[is.na(lmat)] <- 0
    if(verbose) {
	cat("layout: widths = ", lwid, ", heights = ", lhei,"; lmat=\n")
	print(lmat)
    }

    ## Graphics `output' -----------------------

    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
    ## draw the side bars
    if(!missing(RowSideColors)) {
	par(mar = c(margins[1],0, 0,0.5))
	image(rbind(1L:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if(!missing(ColSideColors)) {
	par(mar = c(0.5,0, 0,margins[2]))
	image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    ## draw the main carpet
    par(mar = c(margins[1], 0, 0, margins[2]))
    if(!symm)
	x <- t(x)
    if(revC) { # x columns reversed
	iy <- nr:1
	ddr <- rev(ddr)
	x <- x[,iy]
    } else iy <- 1L:nr

    image(1L:nc, 1L:nr, x, xlim = 0.5+ c(0, nc), ylim = 0.5+ c(0, nr),
	  axes = FALSE, xlab = "", ylab = "", zlim=c(0,2*pi), ...)
    xna <- is.na(x)
    mode(xna) <- "numeric"
    xna[xna==0] <- NA
    image(1L:nc, 1L:nr, xna, col=NAColors, add=TRUE)
    axis(1, 1L:nc, labels= labCol, las= 2, line= -0.5, tick= 0, cex.axis= cexCol)
    if(!is.null(xlab)) mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels= labRow, las= 2, line= -0.5, tick= 0, cex.axis= cexRow)
    if(!is.null(ylab)) mtext(ylab, side = 4, line = margins[2] - 1.25)

    if (!missing(add.expr))
	eval(substitute(add.expr))

    ## the two dendrograms :
    par(mar = c(margins[1], 0, 0, 0))
    if(doRdend)
	plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()

    par(mar = c(0, 0, if(!is.null(main)) 1 else 0, margins[2]))
    if(doCdend)
	plot(ddc,		axes = FALSE, xaxs = "i", leaflab = "none")
    else if(!is.null(main)) frame()
    ## title
    if(!is.null(main)) title(main, cex.main = 1.5*op[["cex.main"]])

    par(mar = c(0, 0, 0, margins[2]))
    if (!is.na(annotate[1]))
	eval(annotate.expr[[1]])
    if (!is.na(annotate[3]))
	eval(annotate.expr[[3]])
    par(mar = c(margins[1], 0, 0, 0))
    if (!is.na(annotate[2]))
	eval(annotate.expr[[2]])
    for (i in 4:length(annotate)) {  
      if (!is.na(annotate[i]))
	eval(annotate.expr[[i]])
    }
    invisible(list(rowInd = rowInd, colInd = colInd,
		   Rowv = if(keep.dendro && doRdend) ddr,
		   Colv = if(keep.dendro && doCdend) ddc ))
}

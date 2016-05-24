plotMAhex <- function (MA, array = 1, xlab = "A", ylab = "M",
                       main = colnames(MA)[array],
                       xlim = NULL, ylim = NULL, status = NULL,
                       values, pch, col, cex, nbin=40,
                       zero.weights = FALSE,
                       style = "colorscale", legend = 1.2, lcex = 1,
                       minarea = 0.04, maxarea = 0.8, mincnt = 2,
                       maxcnt = NULL, trans = NULL, inv = NULL,
                       colorcut = NULL,
                       border = NULL, density = NULL, pen = NULL,
                       colramp = function(n){ LinGray(n,beg = 90,end = 15) },
                       newpage = TRUE, type = c("p", "l", "n"),
                       xaxt = c("s", "n"), yaxt = c("s", "n"),
                       verbose = getOption("verbose"))
{
  if (!requireNamespace("marray", quietly = TRUE))
	stop("cannot process objects without package marray")
  if (!requireNamespace("limma", quietly = TRUE))
	stop("cannot process objects without package limma")
  if(is.null(main))main <- ""
  switch(class(MA),marrayRaw={
        x <- marray::maA(MA[,array])
        y <- marray::maM(MA[,array])
        w <- marray::maW(MA[,array])
      },RGList = {
        MA <- limma::MA.RG(MA[, array])
        array <- 1
        x <- MA$A
        y <- MA$M
        w <- MA$w
    }, MAList = {
        x <- as.matrix(MA$A)[, array]
        y <- as.matrix(MA$M)[, array]
        if (is.null(MA$weights))
            w <- NULL
        else
			w <- as.matrix(MA$weights)[, array]
    }, list = {
        if (is.null(MA$A) || is.null(MA$M))
            stop("No data to plot")
        x <- as.matrix(MA$A)[, array]
        y <- as.matrix(MA$M)[, array]
        if (is.null(MA$weights))
            w <- NULL
        else
			w <- as.matrix(MA$weights)[, array]
    }, MArrayLM = {
        x <- MA$Amean
        y <- as.matrix(MA$coefficients)[, array]
        if (is.null(MA$weights))
            w <- NULL
        else
			w <- as.matrix(MA$weights)[, array]
    }, matrix = {
        narrays <- ncol(MA)
        if (narrays < 2)
            stop("Need at least two arrays")
        if (narrays > 5)
            x <- apply(MA, 1, median, na.rm = TRUE)
        else
			x <- rowMeans(MA, na.rm = TRUE)
        y <- MA[, array] - x
        w <- NULL
    }, ExpressionSet = {
		if (!requireNamespace("Biobase", quietly = TRUE))
          stop("cannot process ExpressionSet objects without package Biobase")
        narrays <- ncol(Biobase::exprs(MA))
        if (narrays < 2)
            stop("Need at least two arrays")
        if (narrays > 5)
            x <- apply(Biobase::exprs(MA), 1, median, na.rm = TRUE)
        else
			x <- rowMeans(Biobase::exprs(MA), na.rm = TRUE)
        y <- Biobase::exprs(MA)[, array] - x
        w <- NULL
        if (missing(main))
            main <- colnames(Biobase::exprs(MA))[array]
    }, AffyBatch = {
		if (!requireNamespace("Biobase", quietly = TRUE) ||
			!requireNamespace("affy", quietly = TRUE))
          stop("cannot process AffyBatch objects without package Biobase and affy")
        narrays <- ncol(Biobase::exprs(MA))
        if (narrays < 2)
            stop("Need at least two arrays")
        if (narrays > 5)
            x <- apply(log2(Biobase::exprs(MA)), 1, median, na.rm = TRUE)
        else
			x <- rowMeans(log2(Biobase::exprs(MA)), na.rm = TRUE)
        y <- log2(Biobase::exprs(MA)[, array]) - x
        w <- NULL
        if (missing(main))
            main <- colnames(Biobase::exprs(MA))[array]
    }, stop("MA is invalid object"))
    if (!is.null(w) && !zero.weights) {
        i <- is.na(w) | (w <= 0)
        y[i] <- NA
    }
    if (is.null(xlim))
        xlim <- range(x, na.rm = TRUE)
    if (is.null(ylim))
        ylim <- range(y, na.rm = TRUE)

    hbin <- hexbin(x,y,xbins=nbin,xbnds=xlim,ybnds=ylim, IDs = TRUE)
    hp <- plot(hbin, legend=legend, xlab = xlab, ylab = ylab, main = main,
               type='n', newpage=newpage)
    ## plot the hexagons
    pushHexport(hp$plot.vp)
    if(is.null(maxcnt)) maxcnt <- max(hbin@count)
    if(is.null(colorcut)) colorcut<-seq(0, 1, length = min(17, maxcnt))
    grid.hexagons(hbin, style=style, minarea = minarea, maxarea = maxarea,
             mincnt = mincnt, maxcnt= maxcnt, trans = trans,
             colorcut = colorcut, density = density, border = border,
             pen = pen, colramp = colramp)
    if (is.null(status) || all(is.na(status))) {
        if (missing(pch))
            pch <- 16
        if (missing(cex))
            cex <- 0.3
        if (missing(col)) {
          clrs <- colramp(length(colorcut)-1)
          col <- clrs[1]
        }
        pp <- inout.hex(hbin,mincnt)
        grid.points(x[pp], y[pp], pch = pch[[1]],
                    gp=gpar(cex = cex[1], col=col, fill=col))
    }
    else {
        if (missing(values)) {
            if (is.null(attr(status, "values")))
                values <- names(sort(table(status), decreasing = TRUE))
            else 
				values <- attr(status, "values")
        }
        sel <- !(status %in% values)
        nonhi <- any(sel)
        if (nonhi) grid.points(x[sel], y[sel], pch = 16, gp=gpar(cex = 0.3))
        nvalues <- length(values)
        if (missing(pch)) {
            if (is.null(attr(status, "pch")))
                pch <- rep(16, nvalues)
            else
				pch <- attr(status, "pch")
        }
        if (missing(cex)) {
            if (is.null(attr(status, "cex"))) {
                cex <- rep(1, nvalues)
                if (!nonhi)
                  cex[1] <- 0.3
            }
            else
				cex <- attr(status, "cex")
        }
        if (missing(col)) {
            if (is.null(attr(status, "col"))) {
                col <- nonhi + 1:nvalues
            }
            else
				col <- attr(status, "col")
        }
        pch <- rep(pch, length = nvalues)
        col <- rep(col, length = nvalues)
        cex <- rep(cex, length = nvalues)
        for (i in 1:nvalues) {
            sel <- status == values[i]
            grid.points(x[sel], y[sel], pch = pch[[i]], gp=gpar(cex = cex[i], col = col[i]))
        }
    }
    popViewport()
	if (legend > 0) {
        inner <- getPlt(hp$plot.vp, ret.unit="inches", numeric=TRUE)[1]
        inner <- inner/hbin@xbins
        ysize <- getPlt(hp$plot.vp, ret.unit="inches", numeric=TRUE)[2]
        pushViewport(hp$legend.vp)
        grid.hexlegend(legend, ysize=ysize, lcex = lcex, inner = inner,
                       style= style, minarea= minarea, maxarea= maxarea,
                       mincnt= mincnt, maxcnt= maxcnt,
                       trans=trans, inv=inv,
                       colorcut = colorcut,
                       density = density, border = border, pen = pen,
                       colramp = colramp)

            #if (is.list(pch))
            #    legend(x = xlim[1], y = ylim[2], legend = values,
            #      fill = col, col = col, cex = 0.9)
            #else legend(x = xlim[1], y = ylim[2], legend = values,
            #    pch = pch, , col = col, cex = 0.9)
        popViewport()
      }
    invisible(list(hbin = hbin, plot.vp = hp$plot.vp, legend.vp = hp$legend.vp))
}

hexMA.loess <- function(pMA, span = .4, col = 'red', n = 200)
{
  fit <- hexVP.loess(pMA$hbin, pMA$plot.vp, span = span, col = col, n = n)
  invisible(fit)
}

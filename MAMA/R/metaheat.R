metaheat<-function (x, Rowv=NA, Colv=NA, distfun = dist, hclustfun = hclust, 
     na.rm = TRUE, do.dendro = c(TRUE, TRUE), legend = 0, legfrac = 8, col = heat.colors(12), r.cex = NULL, c.cex = NULL, ...) 
{
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (is.null(r.cex)) r.cex <- -0.1 + 1/log10(nr)
    if (is.null(c.cex)) c.cex <- -0.4 + 1/log10(nc)
    if (!identical(Rowv, NA)) {
        if (!inherits(Rowv, "dendrogram")) {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            ddr <- reorder(ddr, Rowv)
        }
        else ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else {
        rowInd = 1:nr
        do.dendro[1] = FALSE
    }
    if (!identical(Colv, NA)) {
        if (!inherits(Colv, "dendrogram")) {
            hcc <- hclustfun(distfun(t(x)))
            ddc <- as.dendrogram(hcc)
            ddc <- reorder(ddc, Colv)
        }
        else ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else {
        colInd = 1:nc
        do.dendro[2] = FALSE
    }
    x <- x[rowInd, colInd]
    

    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    do.xaxis = !is.null(colnames(x))
    do.yaxis = !is.null(rownames(x))
    margin = rep(0, 4)
    margin[1] = if (do.xaxis) 
        5     else 2
    margin[2] = if (do.dendro[1]) 
        0     else 2
    margin[3] = if (do.dendro[2]) 
        0     else 2
    margin[4] = if (do.yaxis) 
        5     else 2
    if (do.dendro[1] & do.dendro[2]) {
        ll = matrix(c(0, 3, 2, 1), 2, 2, byrow = TRUE)
        ll.width = c(1, 4)
        ll.height = c(1, 4)
    }
    else if (do.dendro[1]) {
        ll = matrix(c(2, 1), 1, 2, byrow = TRUE)
        ll.width = c(1, 4)
        ll.height = 4
    }
    else if (do.dendro[2]) {
        ll = matrix(c(2, 1), 2, 1, byrow = FALSE)
        ll.width = 4
        ll.height = c(1, 4)
    }
    else {
        ll = matrix(1, 1, 1)
        ll.width = 1
        ll.height = 1
    }
    if (legend %in% 1:4) {
        plotnum = max(ll) + 1
        nc = ncol(ll)
        nr = nrow(ll)
        if (legend == 1) {
            ll = rbind(ll, if (nc == 1) 
                plotnum
            else c(0, plotnum))
            ll.height = c(ll.height, sum(ll.height)/(legfrac - 
                1))
            leg.hor = TRUE
        }
        else if (legend == 2) {
            ll = cbind(if (nr == 1) 
                plotnum
            else c(0, plotnum), ll)
            ll.width = c(sum(ll.width)/(legfrac - 1), ll.width)
            leg.hor = FALSE
        }
        else if (legend == 3) {
            ll = rbind(if (nc == 1) 
                plotnum
            else c(0, plotnum), ll)
            ll.height = c(sum(ll.height)/(legfrac - 1), ll.height)
            leg.hor = TRUE
        }
        else if (legend == 4) {
            ll = cbind(ll, if (nr == 1) 
                plotnum 
            else c(0, plotnum))
            ll.width = c(ll.width, sum(ll.width)/(legfrac - 1))
            leg.hor = FALSE
        }
    }
    layout(ll, widths = ll.width, heights = ll.height, respect = TRUE)
    par(mar = margin)
    image(1:ncol(x), 1:nrow(x), t(x), axes = FALSE, xlab = "", 
        ylab = "", col = col, ...)
    if (do.xaxis) {
        axis(1, 1:ncol(x), las = 2, line = -0.5, tick = 0, labels = colnames(x), 
            cex.axis = c.cex)
    }
    if (do.yaxis) {
        axis(4, 1:nrow(x), las = 2, line = -0.5, tick = 0, labels = rownames(x), 
            cex.axis = r.cex)
    }
    if (do.dendro[1]) {
        mm = margin
        mm[2] = 3
        mm[4] = 0
        par(mar = mm)
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    if (do.dendro[2]) {
        mm = margin
        mm[1] = 0
        mm[3] = 3
        par(mar = mm)
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    if (legend %in% 1:4 & all(unique(as.vector(x)) %in% c(0,1))) {
        dum<-t(matrix(c(1,1)))
        colnames(dum)<-c("neg","pos")
        if (leg.hor) {
            par(mar = c(2.2, margin[2], 1.4, margin[4]))
		#par(mar = c(2.4, margin[2], 1.6, margin[4]))
            barplot(as.table(dum), axes=FALSE, col=col, beside=TRUE, border=NA, cex.names=0.9)
        }
        else {
            par(mar = c(margin[1], 1.7, margin[3], 2))
            barplot(as.table(dum), axes=FALSE, col=col, horiz=TRUE, beside=TRUE, border=NA, cex.names=0.9)

        }
	}

if (legend %in% 1:4 & !all(unique(as.vector(x)) %in% c(0,1))) {
	dummy.x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length = length(col))
	dummy.z <- matrix(dummy.x, ncol = 1)
	if (leg.hor) {
	par(mar = c(2, margin[2], 2, margin[4]))
	image(x = dummy.x, y = 1, z = dummy.z, yaxt = "n", 
		col = col, xlim=c(min(x, na.rm = TRUE),max(x, na.rm = TRUE)))
	}
	else {
	par(mar = c(margin[1], 2, margin[3], 2))
	image(x = 1, y = dummy.x, z = t(dummy.z), xaxt = "n", 
		col = col, ylim=c(min(x, na.rm = TRUE),max(x, na.rm = TRUE)))
     }
     }

    invisible(list(rowInd = rowInd, colInd = colInd))
}
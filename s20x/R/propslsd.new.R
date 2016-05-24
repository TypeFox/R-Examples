propslsd.new = function (crosstablist, conf.level = 0.95, arrowlength = 0.1, ...){
    x = crosstablist$Totals[1:nrow(crosstablist$whole.props), 
        1:ncol(crosstablist$whole.props)]
    nrows = nrow(x)
    ncols = ncol(x)
    qval=abs(qnorm((1-conf.level)/2))
    totals = crosstablist$Totals[, ncol(crosstablist$Totals)][1:nrow(crosstablist$whole.props)]
    phat = crosstablist$row.props

    pm = qval * sqrt(sweep(phat * (1 - phat), 1, totals, "/")/2)
    lsd.low = phat - pm
    lsd.up = phat + pm
    maxy.plot = max(lsd.up)
    plus.min = 0.3
    xmax = (1 + plus.min) * ncols * nrows - plus.min
    xmin = 1 - 2 * plus.min
    op = par(mar = c(5.1, 4.1, 4.1, 3.1))
    plot(c(xmin, xmax), c(0, 1.1 * maxy.plot), type = "n", bty = "n", 
        axes = FALSE, xlab = "", ylab = "Proportion", main = "LSD-display intervals", 
        cex.lab = 1)
    axis(2, las = 2)
    axis(4, las = 2)
    xlabels = if (!is.null(rownames(x))) 
        substr(rownames(x), 1, 7)
    else paste("Row", 1:nrows)
    labels = rep(xlabels, ncols)
    a = matrix(1:(nrows * ncols), nrow = nrows)
    for (j in 1:ncols) {
        start = (j - 1) * nrows * (1 + plus.min)
        if (!is.null(colnames(x)[j])) 
            text(start + median(1:nrows), 1.1 * maxy.plot, colnames(x)[j], 
                adj = 0.5)
        else text(start + median(1:nrows), 1.1 * maxy.plot, paste("Col", 
            j), adj = 0.5)
        if (j > 1) 
            lines(c(start, start), c(0, maxy.plot))
        for (i in 1:nrows) {
            rect(start + i - plus.min, 0, start + i + plus.min, 
                phat[i, j], col = i + 1, border = 1)
            a[i, j] = start + i - plus.min + 0.3
            
            if(lsd.up[i, j] != lsd.low[i, j]){
               arrows(start + i, lsd.low[i, j], start + i, lsd.up[i, j], code = 3, length = arrowlength)
            }
        }
    }
    at = as.vector(a)
    axis(1, at = at, labels = labels, las = 2, cex.axis = 0.8)
    mtext(names(dimnames(crosstablist$Totals))[1], side = 1, 
        line = 4)
    mtext("Row Proportions")
    box()
    invisible(list(totals = totals, phat = phat, pm = pm, lsd.low = lsd.low, 
        lsd.up = lsd.up))
}


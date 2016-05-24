polys <-
function (x, scale = TRUE, labels = dimnames(x)[[1]],
    locations = NULL, nrow = NULL, ncol = NULL, key.loc = NULL,
    key.labels = dimnames(x)[[2]], key.xpd = TRUE, xlim = NULL,
    ylim = NULL, flip.labels = NULL, factx=1,facty=1,
    col.stars = NA, axes = FALSE, frame.plot = axes, main = NULL,
    sub = NULL, xlab = "", ylab = "", cex = 0.8, lwd = 1.1,
    lty = par("lty"), xpd = FALSE, mar = pmin(par("mar"), 1.1 +
        c(2 * axes + (xlab != ""), 2 * axes + (ylab != ""), 1,
            0)), add = FALSE, plot = TRUE, ...)
{
# draw polygons as multivariate graphics
#
    if (is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        stop("'x' must be a matrix or a data frame")
    if (!is.numeric(x))
        stop("data in 'x' must be numeric")
    n.loc <- nrow(x)
    n.seg <- ncol(x)
    if (is.null(locations)) {
        if (is.null(nrow))
            nrow <- ceiling(if (!is.numeric(ncol)) sqrt(n.loc) else n.loc/ncol)
        if (is.null(ncol))
            ncol <- ceiling(n.loc/nrow)
        if (nrow * ncol < n.loc)
            stop("nrow * ncol <  number of observations")
        ff <- if (!is.null(labels)) 2.3
        else 2.1
        locations <- expand.grid(ff * 1:ncol, ff * nrow:1)[1:n.loc, ]
        if (!is.null(labels) && (missing(flip.labels) || !is.logical(flip.labels)))
            flip.labels <- ncol * mean(nchar(labels, type = "c")) > 30
    }
    else {
        if (is.numeric(locations) && length(locations) == 2) {
            locations <- cbind(rep.int(locations[1], n.loc),
                rep.int(locations[2], n.loc))
            if (!missing(labels) && n.loc > 1)
                warning("labels do not make sense for a single location")
            else labels <- NULL
        }
        else {
            if (is.data.frame(locations))
                locations <- data.matrix(locations)
            if (!is.matrix(locations) || ncol(locations) != 2)
                stop("'locations' must be a 2-column matrix.")
            if (n.loc != nrow(locations))
                stop("number of rows of 'locations' and 'x' must be equal.")
        }
        if (missing(flip.labels) || !is.logical(flip.labels))
            flip.labels <- FALSE
    }
    xloc <- locations[, 1]
    yloc <- locations[, 2]
    if (scale) {
        x <- apply(x, 2, function(x) (x - min(x, na.rm = TRUE))/diff(range(x,
            na.rm = TRUE)))
    }
    x[is.na(x)] <- 0
    mx <- max(x <- x * factx)
    if (is.null(xlim))
        xlim <- range(xloc) + c(-mx, mx)
    if (is.null(ylim))
        ylim <- range(yloc) + c(-mx, mx)
    op <- par(mar = mar, xpd = xpd)
    on.exit(par(op))
    if (!add)
        plot(0, type = "n", ..., xlim = xlim, ylim = ylim, main = main,
            sub = sub, xlab = xlab, ylab = ylab, asp = 1, axes = axes)
    if (!plot)
        return()
    s.x <- seq(-factx*n.seg/2,factx*n.seg/2, length=n.seg) 

        for (i in 1:n.loc) {
            xco <- s.x+xloc[i]
            yco <- yloc[i]+x[i,]*facty
            lines(xco,yco,lwd=lwd, ...)
            lines(c(xco[1],xco[length(xco)]),c(yloc[i],yloc[i]),lwd=0.7)
        }
    if (!is.null(labels)) {
        y.off <- mx
        if (flip.labels)
            y.off <- y.off + cex * par("cxy")[2] * ((1:n.loc)%%2 - 0.4)
        text(xloc, yloc - y.off, labels, cex = cex, adj = c(0.5, 1))
    }
    if (!is.null(key.loc)) {
        par(xpd = key.xpd)

        xco <- 2.2*s.x + key.loc[1]
        yco <- key.loc[2] +x[3,]*facty*2.2
        lines(xco,yco,lwd=lwd, ...)
        #lines(c(xco[1],xco[length(xco)]),c(key.loc[1],key.loc[2]))
        segments(xco[1],key.loc[2],xco[length(xco)],key.loc[2],lwd=0.7)

        label.x <- 2.2*s.x + key.loc[1]
        label.y <- key.loc[2] - y.off*1.5
        text(label.x, label.y, labels = key.labels, cex = cex)
    }
    if (frame.plot)
        box(...)
    invisible(locations)
}


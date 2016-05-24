##
filled.contour.poly <- function (x = seq(min(poly[,1]), max(poly[,1]), len = nrow(z)),
               y = seq(min(poly[,2]), max(poly[,2]), len = ncol(z)), 
               z, poly, xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE), 
               zlim = range(z, finite = TRUE), 
               levels = pretty(zlim, nlevels), nlevels = 10,
               color.palette = risk.colors, 
               col = color.palette(length(levels) - 1),
               llevels = levels, labels = NULL, labcex = 0.6,
               drawlabel = TRUE, method = "flattest",
               vfont = c("sans serif", "plain"),
               lcol = par("fg"), lty = par("lty"), lwd = par("lwd"),        
               plot.title, plot.axes, key.title, key.axes, asp = NA, 
               xaxs = "i", yaxs = "i", las = 1, axes = TRUE, ...) 
{
    if(!missing(poly)){
        if (missing(z)) {
            if (!missing(x)) {
                if (is.list(x)) {
                    z <- x$z
                    y <- x$y
                    x <- x$x
                } else {
                    z <- x
                    x <- seq(min(poly[,1]), max(poly[,1]), len = nrow(z))
                }
            } else stop("no `z' matrix specified")
        } else if (is.list(x)) {
            y <- x$y
            x <- x$x
        }
    } else { ## missing poly
        if (missing(z)) {
            if (!missing(x)&&!missing(y)) {
                poly <- y
                if (is.list(x)) {
                    z <- x$z
                    y <- x$y
                    x <- x$x
                } else {
                    z <- x
                    x <- seq(min(poly[,1]), max(poly[,1]), len = nrow(z))
                    y <- seq(min(poly[,2]), max(poly[,2]), len = ncol(z))
                }
            } else stop("no `z' and `poly' matrices specified")
        } else if (is.list(x)) {
            poly <- y
            y <- x$y
            x <- x$x
        }          
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
        stop("increasing x and y values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
    par(las = las)
    mar <- mar.orig
    mar[4] <- mar[2]
    mar[2] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
                yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1], col = col)
    if (missing(key.axes)) {
        if (axes) 
            axis(4)
    }
    else key.axes
    box()
    if (!missing(key.title)) 
        key.title
    mar <- mar.orig
    mar[4] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
        stop("no proper `z' matrix specified")
    if (!is.double(z)) 
        storage.mode(z) <- "double"
    #.Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
    #                        col = col))
    .filled.contour(as.double(x), as.double(y), z, levels=as.double(levels), 
       col = col)
    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            axis(1)
            axis(2)
        }
    }
    else plot.axes
    contour(x, y, z, nlevels, levels, labels,
            xlim, ylim, zlim,
            labcex, drawlabel, method,
            vfont, axes = FALSE, frame.plot = FALSE,
            lcol,  lty, lwd, add=T)
    polygon(poly)
    box()
    if (missing(plot.title)) 
        title(...)
    else plot.title
    invisible()
}

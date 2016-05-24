"biv.plot" <- function(dfxy, br = 10, points = TRUE, density = TRUE,
                       kernel = TRUE, o.include = FALSE, pch,
                       cex, col, h, sub,
                       side = c("top", "bottom", "none"), ...)
{
    ## Verifications
    side <- match.arg(side)
    if (!inherits(dfxy, "data.frame"))
        stop("dfxy should be a data frame")
    if (ncol(dfxy) < 2)
        stop("dfxy should have at least two columns")
    if (missing(pch))
        pch <- 16
    if (missing(cex))
        cex <- 0.5
    if (missing(col))
        col <- grey(0.7)

    ## Graphical parameters
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    lay <- layout(matrix(c(2,4,1,3),2,2, byrow = TRUE), c(3,1),
                  c(1,3), TRUE)

    ## preparation of the data
    x <- dfxy[, 1]
    y <- dfxy[, 2]
    xr <- diff(range(x))
    yr <- diff(range(y))
    xby <- xr/(br-1)
    yby <- yr/(br-1)
    xp <- 0
    xn <- 0
    yp <- 0
    yn <- 0
    if (max(x)>0)
        xp <- seq(0, max(x)+xby, by = xby)
    if (max(y)>0)
        yp <- seq(0, max(y)+yby, by = yby)
    if (min(x)<0)
        xn <- seq(0, min(x)-xby, by = -xby)
    if (min(y)<0)
        yn <- seq(0, min(y)-yby, by = -yby)

    xbr <- c(rev(xn[-1]), xp)
    ybr <- c(rev(yn[-1]), yp)

    ## Cuts the points into classes
    xhist <- hist(x, plot = FALSE, br = xbr)
    yhist <- hist(y, plot = FALSE, br = ybr)

    ## Limits of the graphs
    if (o.include) {
        xlim <- c(min(x, 0)-xr*0.05, max(x, 0)+xr*0.05)
        ylim <- c(min(y, 0)-yr*0.05, max(y, 0)+yr*0.05)
    }
    else {
        xlim <- c(min(x)-xr*0.05, max(x)+xr*0.05)
        ylim <- c(min(y)-yr*0.05, max(y)+yr*0.05)
    }
    xhistlim <- c(0, max(xhist$density)*1.05)
    yhistlim <- c(0, max(yhist$density)*1.05)

    ## The main graph
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    ## background
    plot.default(min(x), min(y), type = "n", xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim,
                 xaxs = "i", yaxs = "i", frame.plot = FALSE)
    abline(v = xbr, col = grey(0.9))
    abline(h = ybr, col = grey(0.9))
    abline(h = 0)
    abline(v = 0)

    ## adds the points
    if(points)
        points(x, y, pch = pch, cex = cex)

    ## an eventual 2D smoothing
    if(kernel) {
        if (missing(h))
            h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
        dens <- MASS::kde2d(x, y, h = h, lims = c(xlim, ylim))
        contour(dens, drawlabels = FALSE, col = col, add = TRUE)
    }

    ## the scale box
    if (side != "none") {
        tra <- paste(" dx = ", signif(xby, 2), " ", "\n", " dy = ",
                     signif(yby, 2), " ", sep = "")
        wt <- strwidth(tra, cex = 1)
        ht <- strheight(tra, cex = 1) * 1.5
        xl <- par("usr")[1]
        yu <- par("usr")[4]
        yd <- par("usr")[3]
        if (side == "top") {
            rect(xl, yu - ht, xl + wt, yu, col = "white", border = 0)
            text(xl + wt/2, yu - ht/2, tra, cex = 1)
        }
        if (side == "bottom") {
            rect(xl, yd + ht, xl + wt, yd, col = "white", border = 0)
            text(xl + wt/2, yd + ht/2, tra, cex = 1)
        }
    }
    box()

    ## Histogram and density of x
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    if(density) {
        xdens <- density(x)
        xhistlim <- c(0, max(xhist$density, xdens$y)*1.05)
    }
    plot.default(min(x), 0, type = "n", xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n", xlim = xlim, ylim = xhistlim,
                 xaxs = "i", yaxs = "i", frame.plot = FALSE)
    rect(xbr[-length(xbr)], rep(0, br), xbr[-1], xhist$density)
    if(density)
        lines(xdens, col = col)
    abline(h = 0)

    ## Histogram and density of y
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    if(density) {
        ydens <- density(y)
        yhistlim <- c(0, max(yhist$density, ydens$y)*1.05)
    }
    plot.default(min(x), 0, type = "n", xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n", xlim = yhistlim, ylim = ylim,
                 xaxs = "i", yaxs = "i", frame.plot = FALSE)
    rect(rep(0, br), ybr[-length(ybr)], yhist$density, ybr[-1])
    if(density)
        lines(ydens$y, ydens$x, col = col)


    ## Main title of the graph
    plot.default(0, 0, type = "n", xlab = "", ylab = "", xaxt = "n",
                 yaxt = "n", xaxs = "i", yaxs = "i", frame.plot = FALSE)
    if (missing(sub))
        sub <- "Biplot and \nmarginals\ndistributions"
    mtext(text = paste(sub), adj = 0.5, line = -8, cex = 1.5)
}

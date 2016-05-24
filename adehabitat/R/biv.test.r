"biv.test" <- function(dfxy, point, br = 10, points = TRUE, density = TRUE,
                       kernel = TRUE, o.include = FALSE, pch, cex,
                       col, Pcol, h, sub,
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
        col <- grey(0.6)
    if (missing(Pcol))
        Pcol <- grey(0.6)

    ## Graphical settings
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    lay <- layout(matrix(c(2,4,1,3),2,2, byrow = TRUE), c(3,1),
                  c(1,3), TRUE)
    layout.show(lay)

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
    if (max(x)>0 | point[1]>0)
        xp <- seq(0, max(x, point[1])+xby, by = xby)
    if (max(y)>0 | point[2]>0)
        yp <- seq(0, max(y, point[2])+yby, by = yby)
    if (min(x)<0 | point[1]<0)
        xn <- seq(0, min(x, point[1])-xby, by = -xby)
    if (min(y)<0 | point[2]<0)
        yn <- seq(0, min(y, point[2])-yby, by = -yby)
    xbr <- c(rev(xn[-1]), xp)
    ybr <- c(rev(yn[-1]), yp)

    ## Cuts the points into classes
    xhist <- hist(x, plot = FALSE, br = xbr, freq = FALSE)
    yhist <- hist(y, plot = FALSE, br = ybr, freq = FALSE)

    ## Limits of the graphs
    if (o.include) {
        xlim <- c(min(x, 0, point[1])-xr*0.05, max(x, 0,
                                                   point[1])+xr*0.05)
        ylim <- c(min(y, 0, point[2])-yr*0.05, max(y, 0,
                                                   point[2])+yr*0.05)
    }
    else {
        xlim <- c(min(x, point[1])-xr*0.05, max(x,
                                                point[1])+xr*0.05)
        ylim <- c(min(y, point[2])-yr*0.05, max(y,
                                                point[2])+yr*0.05)
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

    ## lines showing the abscissa of the observation
    lines(c(point[1], xlim[2]), rep(point[2], 2), lty = 3)
    ## lines showing the ordinate of the observation
    lines(rep(point[1], 2), c(point[2], ylim[2]), lty = 3)

    ## The scale box
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
    ## the observation of the test
    points(point[1], point[2], pch = 18, cex = cex*4, col = Pcol)
    box()


    ## Marginal distribution for X
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

    ## observation of the test
    lines(rep(point[1], 2), c(0, max(xhist$density*2/3)),
          col = Pcol, lwd = 2)
    points(point[1], max(xhist$density*2/3), pch = 18, cex = 2,
           col = Pcol)

    ## ... and P-value
    pX <- (sum(dfxy[,1] >= point[1]) + 1)/(length(dfxy[,1]) + 1)
    if (pX > 0.5)
        pX <- (sum(dfxy[,1] <= point[1]) + 1)/(length(dfxy[,1]) + 1)
    mtext(text = paste("p =", round(pX, 3)), side = 3, adj = 1,
          line = -1)


    ## Marginal distribution for y
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
    abline(v = 0)

    ## observation of the test
    lines(c(0, max(yhist$density*2/3)), rep(point[2], 2),
          col = Pcol, lwd = 2)
    points(max(yhist$density*2/3), point[2], pch = 18, cex = 2,
           col = Pcol)

    ## ... and P-value of the univariate test
    pY <- (sum(dfxy[,2] >= point[2]) + 1)/(length(dfxy[,2]) + 1)
    if (pY > 0.5)
        pY <- (sum(dfxy[,2] <= point[2]) + 1)/(length(dfxy[,2]) + 1)
    mtext(text = paste("p =", round(pY, 3)), side = 3, adj = 1,
          line = -1, las = 0)


    ## Main title of the graph
    plot.default(0, 0, type = "n", xlab = "", ylab = "", xaxt = "n",
                 yaxt = "n", xaxs = "i", yaxs = "i", frame.plot = FALSE)
    if (missing(sub))
        sub <- "Biplot and\n univariate\ntests"
    mtext(text = paste(sub), adj = 0.5, line = -8, cex = 1.5)
}

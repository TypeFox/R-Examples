scatterniche <- function (x, pr, xax = 1, yax = 2, pts = FALSE, percent = 95,
                          clabel = 1, side = c("top", "bottom", "none"),
                          Adensity, Udensity, Aangle, Uangle, Aborder,
                          Uborder, Acol, Ucol, Alty,
                          Ulty, Abg, Ubg, Ainch, Uinch, ...)
{
    ## Graphical settings
    side <- match.arg(side)
    opar <- par(mar = c(0.1, 0.1, 0.1, 0.1))
    on.exit(par(opar))

    ## Bases for the graphs
    x1 <- x[,xax]
    x1 <- c(x1 - diff(range(x1)/50), x1 + diff(range(x1))/50)
    xlim <- range(x1)
    y1 <- x[, yax]
    y1 <- c(y1 - diff(range(y1)/50), y1 + diff(range(y1))/50)
    ylim <- range(y1)

    ## background graph
    scatterutil.base(dfxy = x[, c(xax, yax)], xax = 1, yax = 2,
                     xlim = xlim, ylim = ylim, grid = TRUE, addaxes = FALSE,
                     cgrid = 1, include.origin = TRUE,
                     origin = c(0, 0), sub = "",
                     csub = 1.25, possub = "bottomleft",
                     pixmap = NULL, contour = NULL,
                     area = NULL, add.plot = FALSE)

    ## If points are desired
    if (pts) {

        ## graphical settings
        if (missing(Acol))
            Acol <- gray(0.8)
        if (missing(Ucol))
            Ucol <- "black"
        if (missing(Abg))
            Abg <- gray(0.8)
        if (missing(Ubg))
            Ubg <- "black"
        if (missing(Ainch))
            Ainch <- 0.03
        if (missing(Uinch))
            Uinch <- Ainch*max(pr)/min(pr[pr>1e-7])

        ## draws the points
        symbols(x[, c(xax, yax)], circles = rep(1, length(pr)),
                fg = Acol, bg = Abg, inches = Ainch, add = TRUE)
        symbols(x[pr > 0, c(xax, yax)], circles = pr[pr > 0],
                fg = Ucol, bg = Ubg,
                inches = Uinch, add = TRUE)
        abline(v = 0)
        abline(h = 0)
    }
    else {
        ## if polygons are desired

        ## Graphical settings
        if (missing(Adensity))
            Adensity <- NULL
        if (missing(Udensity))
            Udensity <- NULL
        if (missing(Aangle))
            Aangle <- 45
        if (missing(Uangle))
            Uangle <- 45
        if (missing(Aborder))
            Aborder <- NULL
        if (missing(Uborder))
            Uborder <- NULL
        if (missing(Acol))
            Acol <- gray(0.95)
        if (missing(Ucol))
            Ucol <- gray(0.6)
        if (missing(Alty))
            Alty <- NULL
        if (missing(Ulty))
            Ulty <- NULL

        ## Convex polygons
        pcff <- function(xy)
        {
            mo <- apply(xy,2,mean)
            dis <- apply(xy, 1, function(x) sum((x-mo)^2))
            xy <- xy[dis < quantile(dis, percent/100),]
            return(xy[chull(xy[,1], xy[,2]),])
        }
        mcpA <- pcff(x[, c(xax, yax)])
        mcpU <- pcff(x[pr>1e-7, c(xax, yax)])
        polygon(mcpA, density = Adensity, angle = Aangle,
                border = Aborder, col = Acol, lty = Alty)
        polygon(mcpU, density = Udensity, angle = Uangle,
                border = Uborder, col = Ucol, lty = Ulty)
        abline(v = 0)
        abline(h = 0)
    }

    ## The legend
    xax <- paste("Axis", xax)
    yax <- paste("Axis", yax)
    if (side != "none") {
        tra <- paste(" xax =", xax, "\n yax =", yax)
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
}

gx.hist <-
function (xx, xlab = deparse(substitute(xx)), ylab = "Number of Observations", 
    log = FALSE, xlim = NULL, main = "", nclass = NULL, colr = NULL, 
    ifnright = TRUE, cex = 0.8, ...) 
{
    temp.x <- remove.na(xx)
    x <- temp.x$x[1:temp.x$n]
    nx <- temp.x$n
    if ((is.null(nclass)) && (nx <= 500)) nclass <- "scott"
    if ((is.null(nclass)) && (nx > 500)) nclass <- "fd"
    xrange <- range(x)
    if (log) {
        logx <- "x"
        x.save <- x
        x <- log10(x)
        if ((!is.null(xlim)) && (xlim[1] <= 0)) 
            xlim[1] <- xrange[1]
    }
    else logx <- ""
    q <- as.vector(quantile(x, c(0.25, 0.75)))
    h <- switch(pmatch(nclass, c("sturges", "Sturges", "scott", 
        "Scott", "fd", "FD"), nomatch = ""), `1` = diff(range(x))/(logb(nx, 
        2) + 1), `2` = diff(range(x))/(logb(nx, 2) + 1), `3` = 3.5 * 
        sqrt(var(x)) * nx^(-1/3), `4` = 3.5 * sqrt(var(x)) * 
        nx^(-1/3), `5` = 2 * (q[2] - q[1]) * nx^(-1/3), `6` = 2 * 
        (q[2] - q[1]) * nx^(-1/3))
    nhb <- ceiling(diff(range(x))/h)
    breaks <- pretty(x, nhb)
    if (log) {
        breaks <- 10^breaks
        x <- x.save
    }
    nb <- length(breaks)
    nk <- nb - 1
    k <- cut(x, breaks, include.lowest = TRUE)
    nink <- as.vector(table(k))
    kmax <- max(nink)
    if (is.null(xlim)) {
        xlim <- c(breaks[1], breaks[nb])
        nnx <- 0
    }
    else {
        dropped <- x[(x > xlim[2]) | (x < xlim[1])]
        nnx <- length(dropped)
    }
    nvec <- (nk * 2) + 2
    xpos <- numeric(nvec + 1)
    ypos <- numeric(nvec + 1)
    for (i in 1:nb) {
        for (j in 1:2) {
            xpos[(i - 1) * 2 + j] <- breaks[i]
        }
    }
    xpos[nvec + 1] <- breaks[1]
    ypos[1] <- 0
    for (i in 1:nk) {
        for (j in 1:2) {
            ypos[(i - 1) * 2 + j + 1] <- nink[i]
        }
    }
    ypos[nvec] <- 0
    ypos[nvec + 1] <- 0
    if (!is.null(xlim)) {
        x <- xpos[(xpos >= xlim[1]) & (xpos <= xlim[2])]
        y <- ypos[(xpos >= xlim[1]) & (xpos <= xlim[2])]
        nvec <- length(x)
        x[1:2] <- xlim[1]
        x[nvec - 0:1] <- xlim[2]
        x[nvec + 1] <- xlim[1]
        y[1] <- 0
        y[nvec] <- 0
        y[nvec + 1] <- 0
    }
    else {
        x <- xpos
        y <- ypos
    }
    plot(x, y, log = logx, xlim = xlim, ylim = c(0, kmax), xlab = xlab, 
        ylab = ylab, main = main, type = "n", las = 1, ...)
    if(is.null(colr)) colr <- 8
    polygon(x, y, col = colr)
    lines(x, y)
    limits <- par("usr")
    if (!is.null(ifnright)) {
        if (ifnright) {
            xpos <- limits[2] - (limits[2] - limits[1]) * 0.05
            adj <- 1
        }
        else {
            xpos <- limits[1] + (limits[2] - limits[1]) * 0.05
            adj <- 0
        }
        if (log) 
            xpos <- 10^xpos
            ypos <- limits[4] - (limits[4] - limits[3]) * 0.1
            text(xpos, ypos, labels = paste("N =", nx), adj = adj, cex = cex, 
                ...)
    }
    if (nnx >= 1) {
        ypos <- limits[4] - (limits[4] - limits[3]) * 0.2
        text(xpos, ypos, labels = paste("Bins for", nnx, "\npoints omitted"), 
            cex = cex * 0.8, adj = adj, ...)
    }
    invisible(list(xlim = xlim))
}

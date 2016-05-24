caplot <-
function (x, y, z, zname = deparse(substitute(z)), caname = NULL, 
    log = TRUE, ifjit = FALSE, ifrev = FALSE, ngrid = 100, colr = topo.colors(16), 
    xcoord = "Easting", ycoord = "Northing") 
{
    frame()
    oldpar <- par()
    on.exit(par(oldpar))
    par(mfrow = c(2, 2), pty = "m", cex.main = 0.8)
    u <- na.exclude(cbind(x, y, z))
    xlim <- range(u[, 3])
    cnpplt(u[, 3], xlab = zname, log = log, xlim = xlim, main = "% Cumulative Probability Plot\nOriginal Data", 
        cex.axis = 1, ifshape = TRUE, cex.lab = 0.8)
    if (ifjit) {
        u[, 1] <- jitter(u[, 1], 0.5)
        u[, 2] <- jitter(u[, 2], 0.5)
    }
        if (log) u[, 3] <- log10(u[, 3])
        if (is.null(caname)) {
            if (log) caname <- paste("Log10(", zname, ")", sep = "")
            else caname <- zname 
    }
    xo <- seq(min(u[, 1]), max(u[, 1]), length.out = ngrid)
    yo <- seq(min(u[, 2]), max(u[, 2]), length.out = ngrid)
    new <- akima::interp(u[, 1], u[, 2], u[, 3], xo = xo, yo = yo)
    znew <- na.omit(as.vector(new$z))
    if (log) 
        znew <- 10^znew
    cnpplt(as.vector(znew), xlab = zname, log = log, xlim = xlim, 
        main = "% Cumulative Probability Plot\nGridded Data", 
        cex.axis = 0.8, ifshape = TRUE, cex.lab = 0.8)
    eqscplot(range(new$x), range(new$y), plot = "n", xlab = xcoord, 
        ylab = ycoord, main = caname, pch = 32, cex.lab = 0.8)
    image(new, add = TRUE, col = colr)
    gx.mf(znew, ifrev, xlab = zname, ylab = "Cumulative Area (%)",
        main = "Concentration-Area (C-A) Plot", xlim = xlim, 
        cex.lab = 0.8)
    invisible()
}

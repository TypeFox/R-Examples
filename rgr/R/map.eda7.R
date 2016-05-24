map.eda7 <-
function (xx, yy, zz, sfact = 1, logz = FALSE, xlab = "Easting", 
    ylab = "Northing", zlab = deparse(substitute(zz)), main = "", 
    ifgrey = FALSE, symcolr = NULL, tol = 0.04, iflgnd = FALSE, 
    title = deparse(substitute(zz)), cex.lgnd = 0.8, ...) 
{
    frame()
    oldpar <- par()
    on.exit(par(oldpar))
    par(pty = "m")
    temp.x <- remove.na(cbind(xx, yy, zz))
    x <- temp.x$x[1:temp.x$n, 1]
    y <- temp.x$x[1:temp.x$n, 2]
    z <- temp.x$x[1:temp.x$n, 3]
    nz <- temp.x$n
    if (main == "") 
        if (zlab == "") 
            banner <- ""
        else {
            if (logz) 
                banner <- paste("EDA Tukey Boxplot Based Map for Log10", 
                  zlab)
            else banner <- paste("EDA Tukey Boxplot Based Map for", 
                zlab)
        }
    else banner <- main
    eqscplot(x, y, type = "n", xlab = xlab, ylab = ylab, main = banner, 
        tol = tol, ...)
    if (logz) 
        z <- log10(z)
    q <- quantile(z, probs = c(0.25, 0.75))
    zcut <- numeric(6)
    hw <- q[2] - q[1]
    zcut[1] <- q[1] - 3 * hw
    zcut[2] <- q[1] - 1.5 * hw
    zcut[3] <- q[1]
    zcut[4] <- q[2]
    zcut[5] <- q[2] + 1.5 * hw
    zcut[6] <- q[2] + 3 * hw
    zzz <- cutter(z, zcut)
    npch <- c(1, 1, 1, 3, 0, 0, 0)
    size <- c(2, 1.5, 1, 0.5, 1, 1.5, 2) * sfact
    if (ifgrey) {
        symcolr <- grey(c(0, 0.15, 0.3, 0.4, 0.3, 0.15, 0))
    }
    else {
        palette(rainbow(36))
        if (length(symcolr) != 7) 
            symcolr <- c(25, 22, 20, 13, 6, 4, 1)
    }
    for (i in 1:nz) {
        points(x[i], y[i], pch = npch[zzz[i]], cex = size[zzz[i]], 
            col = symcolr[zzz[i]])
    }
    cat("\tCut Levels\t  No. of Symbols   Symbol - size - Colour\n\tLog =", 
        logz, "\t\t\t\tsfact =", format(sfact, nsmall = 2), "\n\n")
    stype <- character(7)
    stype[1:3] <- "Circle "
    stype[4] <- "Cross  "
    stype[5:7] <- "Square "
    for (i in 1:6) {
        if (logz) 
            zcut[i] <- 10^zcut[i]
        ni <- length(zzz[zzz == i])
        cat("\t\t\t      ", ni, "\t    ", stype[i], format(size[i], 
            nsmall = 2), "  ", symcolr[i], "\n\t", round(zcut[i], 
            2), "\n")
    }
    cat("\t\t\t      ", length(zzz[zzz == 7]), "\t    ", stype[7], 
        format(size[7], nsmall = 2), "  ", symcolr[7], "\n")
    if (iflgnd) {
        lgnd.line <- numeric(7)
        zcut <- signif(zcut, 3)
        lgnd.line[1] <- paste(">", zcut[6])
        lgnd.line[2] <- paste(zcut[5], "-", zcut[6])
        lgnd.line[3] <- paste(zcut[4], "-", zcut[5])
        lgnd.line[4] <- paste(zcut[3], "-", zcut[4])
        lgnd.line[5] <- paste(zcut[2], "-", zcut[3])
        lgnd.line[6] <- paste(zcut[1], "-", zcut[2])
        lgnd.line[7] <- paste("<", zcut[1])
        legend(locator(1), pch = npch[7:1], col = symcolr[7:1], 
            pt.cex = size[7:1], lgnd.line[1:7], cex = cex.lgnd,
            title = title, ...)
    }
    palette("default")
    invisible()
}

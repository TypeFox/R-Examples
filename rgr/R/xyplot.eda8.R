xyplot.eda8 <-
function (xx, yy, zz, sfact = 1, xlim = NULL, ylim = NULL, xlab = deparse(substitute(xx)), 
    ylab = deparse(substitute(yy)), zlab = deparse(substitute(zz)), 
    main = "", log = NULL, ifgrey = FALSE, symcolr = NULL, iflgnd = FALSE, 
    pctile = FALSE, title = deparse(substitute(zz)), cex.lgnd = 0.8, ...) 
{
    frame()
    oldpar <- par()
    on.exit(par(oldpar))
    temp.z <- remove.na(cbind(xx, yy, zz))
    x <- temp.z$x[1:temp.z$n, 1]
    y <- temp.z$x[1:temp.z$n, 2]
    z <- temp.z$x[1:temp.z$n, 3]
    nz <- temp.z$n
    if (main == "") 
        if (zlab == "") 
            banner <- ""
        else banner <- paste("EDA Percentile Based Plot for", 
            zlab)
    else banner <- main
    if (is.null(log)) 
        log = ""
    plot(x, y, type = "n", log = log, xlab = xlab, ylab = ylab, 
        xlim = xlim, ylim = ylim, main = banner, ...)
    zcut <- quantile(z, probs = c(0.02, 0.05, 0.25, 0.5, 0.75, 
        0.95, 0.98))
    zzz <- cutter(z, zcut)
    npch <- c(1, 1, 1, 1, 0, 0, 0, 0)
    size <- c(2, 1.5, 1, 0.5, 0.5, 1, 1.5, 2) * sfact
    if (ifgrey) {
        symcolr <- grey(c(0, 0.15, 0.3, 0.4, 0.4, 0.3, 0.15, 
            0))
    }
    else {
        palette(rainbow(36))
        if (length(symcolr) != 8) 
            symcolr <- c(25, 22, 20, 13, 13, 6, 4, 1)
    }
    for (i in 1:nz) {
        points(x[i], y[i], pch = npch[zzz[i]], cex = size[zzz[i]], 
            col = symcolr[zzz[i]])
    }
    cat("\tCut Levels\t No. of Symbols   Symbol - size - Colour\n\t\t\t\t\t\tsfact =", 
        format(sfact, nsmall = 2), "\n\n")
    stype <- character(8)
    stype[1:4] <- "Circle"
    stype[5:8] <- "Square"
    pct <- 0
    for (i in 1:7) {
        ni <- length(zzz[zzz == i])
        pct <- pct + 100 * ni/nz
        cat("\t\t\t      ", ni, "\t    ", stype[i], format(size[i], 
            nsmall = 2), "  ", symcolr[i], "\n\t", signif(zcut[i], 
            4), "\t", round(pct, 1), "%\n")
    }
    ni <- length(zzz[zzz == 8])
    cat("\t\t\t      ", ni, "\t    ", stype[8], format(size[8], 
        nsmall = 2), "  ", symcolr[8], "\n")
    if (iflgnd) {
        lgnd.line <- numeric(8)
        zcut <- signif(zcut, 3)
        if (pctile) {
            title <- paste(deparse(substitute(zz)), "Percentiles")
            lgnd.line[1] <- "> 98th"
            lgnd.line[2] <- "95th - 98th"
            lgnd.line[3] <- "75th - 95th"
            lgnd.line[4] <- "50th - 75th"
            lgnd.line[5] <- "25th - 50th"
            lgnd.line[6] <- "5th - 25th"
            lgnd.line[7] <- "2nd - 5th"
            lgnd.line[8] <- "< 2nd"
        }
        else {
            lgnd.line[1] <- paste(">", zcut[7])
            lgnd.line[2] <- paste(zcut[6], "-", zcut[7])
            lgnd.line[3] <- paste(zcut[5], "-", zcut[6])
            lgnd.line[4] <- paste(zcut[4], "-", zcut[5])
            lgnd.line[5] <- paste(zcut[3], "-", zcut[4])
            lgnd.line[6] <- paste(zcut[2], "-", zcut[3])
            lgnd.line[7] <- paste(zcut[1], "-", zcut[2])
            lgnd.line[8] <- paste("<", zcut[1])
        }
        legend(locator(1), pch = npch[8:1], col = symcolr[8:1], 
            pt.cex = size[8:1], lgnd.line[1:8], cex = cex.lgnd,
            title = title, ...)
    }
    palette("default")
    invisible()
}

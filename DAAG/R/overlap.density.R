overlap.density <-
function (x0, x1, ratio = c(0.05, 20), ratio.number = FALSE, 
          plotvalues = c("Density", "Numbers"), gpnames = c("Control", "Treatment"), 
          cutoffs = c(lower = TRUE, upper = TRUE), bw = FALSE, 
          xlab = "Score", ylab = NULL, 
          col = 1:2, lty = 1:2, lwd = c(1, 1)) 
{
    if (is.null(plotvalues)) 
        plotvalues <- ""
    if (length(plotvalues) > 1) 
        plotvalues <- plotvalues[1]
    if (plotvalues %in% c("Numbers", "Density")) 
        plotit <- TRUE else plotit <- FALSE
    ran <- range(c(x0, x1))
    if (all(cutoffs)) {
        d0 <- density(x0, from = ran[1], to = ran[2])
        d1 <- density(x1, from = ran[1], to = ran[2])
    }
    else if (cutoffs[1]) {
        d0 <- density(x0, from = ran[1])
        d1 <- density(x1, from = ran[1])
    }
    else if (cutoffs[2]) {
        d0 <- density(x0, to = ran[2])
        d1 <- density(x1, to = ran[2])
    }
    else {
        d0 <- density(x0)
        d1 <- density(x1)
    }
    n0 <- length(x0)
    n1 <- length(x1)
    f0 <- d0$y
    f1 <- d1$y
    if (plotvalues %in% "Number") {
        pf0 <- d0$y * n0
        pf1 <- d1$y * n1
    }
    else {
        pf0 <- d0$y
        pf1 <- d1$y
    }
    if (plotit) {
        xlim <- range(c(d0$x, d1$x), na.rm = TRUE)
        if (plotvalues == "Number" & is.null(ylab)) 
            ylab <- "Density x total frequency"
        else if (plotvalues == "Density" & is.null(ylab)) 
            ylab <- "Density"
        ylim <- range(c(0, pf0, pf1))
        ylim[2] <- ylim[2] + 0.1 * diff(ylim)
        plot(d1$x, pf1, xlim = xlim, xlab = xlab, xaxt = "n", 
            bty = "n", yaxs = "i", ylim = ylim, ylab = ylab, 
            main = "", type = "n")
        axis(1)
        box(bty = "L")
        lines(d0$x, pf0, col = col[1], lty = lty[1], lwd = lwd[1])
        if (bw & lty[2] > 1) 
            lines(d1$x, pf1, col = col[2], lty = 1)
        lines(d1$x, f1, col = col[2], lty = lty[2], lwd = lwd[2])
        xpos <- par()$usr[2]
        ypos <- par()$usr[4]
        legend(xpos, ypos, lty = lty, lwd = lwd, col = col, cex = 0.8, 
            legend = gpnames, bty = "n", xjust = 1)
    }
    if (is.null(ratio)) 
        return()
    d0 <- density(x0, from = xlim[1], to = xlim[2])
    d1 <- density(x1, from = xlim[1], to = xlim[2])
    x01 <- d0$x
    f0 <- d0$y
    f1 <- d1$y
    if (ratio.number) {
        cf0 <- f0 * n0
        cf1 <- f1 * n1
    }
    else {
        cf0 <- f0
        cf1 <- f1
    }
    cf0[cf0 < 0] <- 0
    cf1[cf1 < 0] <- 0
    fmin <- pmin(cf0, cf1)
    fmax <- max(fmin)
    subs <- match(fmax, fmin)
    xmid <- x01[subs]
    flow <- ratio[1]
    fhi <- ratio[2]
    lochoose <- x01 < xmid & (cf0 <= flow * cf1 | cf1 <= cf0 * 
        flow)
    if (any(lochoose)) 
        xlim[1] <- max(x01[lochoose])
    else xlim[1] <- min(x01)
    hichoose <- x01 > xmid & (cf0 >= fhi * cf1 | cf1 >= cf0 * 
        fhi)
    if (any(hichoose)) 
        xlim[2] <- min(x01[hichoose])
    else xlim[2] <- max(x01)
    if (plotit) {
        if(ratio.number)midlab <- "ratio: number" else midlab <- "ratio: p.d."
        axis(3, at = xlim, labels = paste(signif(xlim, 4)), mgp = c(3, 
            0.5, 0), col = "gray45", line = 1, cex.axis = 0.8)
        xlim3 <- c(xlim[1], mean(xlim), xlim[2])
        axis(3, at = xlim3, line = 1, labels = c(paste("1:", round(1/ratio[1]), 
            sep = ""), midlab, paste(round(ratio[2]), ":1", sep = "")), 
            mgp = c(-3, -1, 0), col = "gray", cex.axis = 0.7, 
            lty = 0)
    }
    xlim
}

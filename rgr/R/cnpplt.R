cnpplt <-
function (xx, xlab = deparse(substitute(xx)), ylab = "% Cumulative Probability", 
    log = FALSE, xlim = NULL, main = "", ifqs = FALSE, ifshape = FALSE, 
    pch = 3, cex = 0.8, cexp = 1, cex.axis = 0.8, ...) 
{
    temp.x <- remove.na(xx)
    x <- sort(temp.x$x[1:temp.x$n])
    nx <- temp.x$n
    y <- qnorm(ppoints(nx))
    if (log) {
        logx <- "x"
        if ((!is.null(xlim)) && (xlim[1] <= 0)) 
            xlim[1] <- min(x)
    }
    else logx <- ""
    if (is.null(xlim)) {
        plot(x, y, xlab = xlab, ylab = ylab, log = logx, yaxt = "n", 
            main = main, type = "n", ...)
        points(x, y, pch = pch, cex = cexp)
        limits <- par("usr")
        nxx <- nx
    }
    else {
        xt <- x[(x >= xlim[1]) & (x <= xlim[2])]
        yt <- y[(x >= xlim[1]) & (x <= xlim[2])]
        plot(xt, yt, xlab = xlab, ylab = ylab, log = logx, xlim = xlim, 
            yaxt = "n", main = main, type = "n", ...)
        points(xt, yt, pch = pch, cex = cexp)
        limits <- par("usr")
        nxx <- length(xt)
    }
    labels = c("0.001", "0.01", "0.1", "1", "2", "5", "10", "25", 
        "50", "75", "90", "95", "98", "99", "99.9", "99.99", 
        "99.999")
    if (ifshape) {
        if (nxx > 100) {
            labels[1:3] <- ""
            labels[5:6] = ""
            labels[8] <- ""
            labels[10] <- ""
            labels[12:13] = ""
            labels[15:17] <- ""
        }
        else {
            labels[1:3] <- ""
            labels[5] <- ""
            labels[7] <- ""
            labels[11] <- ""
            labels[13] <- ""
            labels[15:17] <- ""
        }
        cex.axis = 0.8
    }
    axis(2, at = c(-4.264891, -3.719016, -3.090232, -2.326348, 
        -2.053749, -1.644854, -1.281552, -0.67449, 0, 0.67449, 
        1.281552, 1.644854, 2.053749, 2.326348, 3.090232, 3.719016, 
        4.264891), labels = labels, las = 1, cex.axis = cex.axis)
    xpos <- limits[2] - (limits[2] - limits[1]) * 0.05
    ypos <- limits[3] + (limits[4] - limits[3]) * 0.11
    if (log) 
        xpos <- 10^xpos
    text(xpos, ypos, labels = paste("N =", nx), adj = 1, cex = cex)
    if (nxx != nx) {
        ypos <- limits[3] + (limits[4] - limits[3]) * 0.05
        text(xpos, ypos, labels = paste(nx - nxx, "points omitted"), 
            adj = 1, cex = cex * 0.8)
    }
    if (ifqs) {
        abline(h = qnorm(c(0.25, 0.5, 0.75)), lty = 3)
        abline(v = quantile(x, probs = c(0.25, 0.5, 0.75)), lty = 3)
    }
    invisible()
}

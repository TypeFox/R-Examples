plot.summary.rmas <-
function (x, ylim = NULL, col = NULL, xlab = NULL, ylab = NULL, 
    main = NULL, ...) 
{
    cosa.s <- x
    if (is.null(ylab)) 
        ylab <- "Abundance"
    if (is.null(xlab)) 
        xlab <- "Time"
    if (is.null(ylim)) 
        ylim <- range(cosa.s[, -1],na.rm=TRUE)
    if (is.null(col)) 
        col <- c("blue", "grey", "red", "red")
    if (dim(x)[2] == 2) {
        plot(cosa.s$Time, cosa.s$Abundance, type = "l", ylim = ylim, 
            col = col[1], ylab = ylab, xlab = xlab, ...)
    }
    if (dim(x)[2] == 6) {
        plot(cosa.s$Time, cosa.s$Average, type = "l", ylim = ylim, 
            col = col[1], ylab = ylab, xlab = xlab, ...)
        arrows(x0 = cosa.s[-1, 1], y0 = cosa.s[-1, 3], x1 = cosa.s[-1, 
            1], y1 = cosa.s[-1, 5], code = 3, angle = 90, length = 0.025, 
            col = col[2])
        points(cosa.s[-1, 1], cosa.s[-1, 2], col = col[3])
        points(cosa.s[-1, 1], cosa.s[-1, 6], col = col[4])
    }
}

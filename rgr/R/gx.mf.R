gx.mf <-
function (xx, xlab = deparse(substitute(xx)), 
    ylab = "Cumulative Percentage of Data", main = "Multifractality Plot",
    ifrev = FALSE, xlim = range(xx, na.rm = TRUE), ...) 
{
    x <- na.exclude(xx)
    nx <- length(x)
    x <- x[order(x)]
    if (!ifrev) x <- rev(x)
    cumn <- seq(1, nx)/nx * 100
    plot(x, cumn, log = "xy", xlab = xlab, ylab = ylab, main = main,
        xlim = xlim, yaxt = "n", pch = 3, ...)
    axis(2, at = c(0.01, 0.1, 1, 10, 100), labels = c("0.01", 
        "0.1", "1", "10", "100"), las = 1, ...)
    limits <- par("usr")
    xpos <- 10^(limits[2] - (limits[2] - limits[1]) * 0.05)
    if (ifrev) 
        ypos <- 10^(limits[3] + (limits[4] - limits[3]) * 0.11)
    else ypos <- 10^(limits[4] - (limits[4] - limits[3]) * 0.11)
    text(xpos, ypos, labels = paste("N =", nx), adj = 1, cex = 0.8)
    invisible()
}

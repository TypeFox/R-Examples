bwplots.by.var <-
function (xmat, log = FALSE, wend = 0.05, notch = FALSE, xlab = "Measured Variables", 
    ylab = "Reported Values", main = "", label = NULL, plot.order = NULL, 
    xpos = NA, las = 1, cex.axis = 1, adj = 0.5, colr = 8, pch = 3, 
    ...) 
{
    zz <- var2fact(xmat)
    x <- zz[, 1]
    y <- as.numeric(zz[, 2])
    if (is.null(label)) 
        label <- sort(unique(x))
    bwplots(split(y, x), log = log, wend = wend, notch = notch, 
        xlab = xlab, ylab = ylab, main = main, label = label, 
        plot.order = plot.order, xpos = xpos, las = las, cex.axis = cex.axis, 
        adj = adj, colr = colr, pch = pch, ...)
    invisible()
}

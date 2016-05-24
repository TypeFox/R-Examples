tbplots.by.var <-
function (xmat, log = FALSE, logx = FALSE, notch = FALSE, xlab = "Measured Variables", 
    ylab = "Reported Values", main = "", label = NULL, plot.order = NULL, 
    xpos = NA, las = 1, cex = 1, adj = 0.5, colr = 8, ...) 
{
    zz <- var2fact(xmat)
    x <- zz[, 1]
    y <- as.numeric(zz[, 2])
    if (is.null(label)) 
        label <- sort(unique(x))
    tbplots(split(y, x), log = log, logx = logx, notch = notch, 
        xlab = xlab, ylab = ylab, main = main, label = label, 
        plot.order = plot.order, xpos = xpos, las = las, cex = cex, 
        adj = adj, colr = colr, ...)
    invisible()
}

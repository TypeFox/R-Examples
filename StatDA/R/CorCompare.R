"CorCompare" <-
function (cor1, cor2, labels1, labels2, method1, method2, ndigits = 4, 
    lty1=1, lty2=2, col1=1, col2=2, lwd1=1.1, lwd2=1.1, 
    cex.label=1.1, cex.legend=1.2, lwd.legend=1.2, cex.cor=1, ...) 
{
# Compare two correlation matrices numerically and graphically
#
# cor1, cor2 ... two correlation matrices based on different estimation methods
# labels1, labels2 ... lables  for the two estimation methods
# method1, method2 ... description of the estimation methods
# rest ... other graphics parameters
#
    plot.new()
    par(mar = rep(1, 4))
    p = ncol(cor1)
    lim = c(-1, p + 1)
    plot.window(xlim = lim, ylim = lim, xaxs = "i", yaxs = "i")
    for (i in 1:p) {
        text(i - 0.5, p + 0.5, labels1[i], srt = 90, cex=cex.label)
        text(-0.5, p - i + 0.5, labels2[i], cex=cex.label)
    }
    for (i in 2:p) {
        for (j in 1:(i - 1)) {
            do.ellipses(cor1[c(i, j), c(i, j)], pos = c(i - 1.5, 
                p - j - 0.5), lwd = lwd1, lty=lty1)
            do.ellipses(cor2[c(i, j), c(i, j)], pos = c(i - 1.5, 
                p - j - 0.5), lwd = lwd2, lty=lty2, ...)
            text(j - 0.5, p - i + 0.5, round(cor1[i, j], ndigits), 
                adj = c(0.5, -0.1), cex=cex.cor)
            text(j - 0.5, p - i + 0.5, round(cor2[i, j], ndigits), 
                adj = c(0.5, 1.1),, cex=cex.cor, ...)
        }
    }
    lines(c(0.5, p - 0.5), c(p - 0.5, 0.5), lwd = 2, lty=3)
    lines(lim[2] - 2 + c(-0.5, 0.3), c(-0.3, -0.3), lwd=lwd.legend, lty=lty1)
    lines(lim[2] - 2 + c(-0.5, 0.3), c(-0.7, -0.7), lwd=lwd.legend, lty=lty2, ...)
    text(lim[2] - 2 - 0.7, -0.275, method1, pos = 2, cex=cex.legend)
    text(lim[2] - 2 - 0.7, -0.675, method2, pos = 2, cex=cex.legend)
}

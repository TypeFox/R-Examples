gx.ks.test <-
function (xx1, xx2, xlab = " ", x1lab = deparse(substitute(xx1)), 
    x2lab = deparse(substitute(xx2)), ylab = "Empirical Cumulative Distribution Function", 
    log = FALSE, main = "", pch1 = 3, col1 = 2, pch2 = 4, col2 = 4, 
    ifresult = TRUE, cex = 0.8, cexp = 0.9, ...) 
{
    temp.x <- remove.na(xx1)
    x1 <- sort(temp.x$x[1:temp.x$n])
    nx1 <- temp.x$n
    y1 <- ((1:nx1) - 0.5)/nx1
    temp.x <- remove.na(xx2)
    x2 <- sort(temp.x$x[1:temp.x$n])
    nx2 <- temp.x$n
    y2 <- ((1:nx2) - 0.5)/nx2
    xlim <- range(c(x1, x2))
    if (log) {
        logx <- "x"
        if (xlim[1] <= 0) 
            stop("\n  Values cannot be .le. zero for a log plot\n")
    }
    else logx <- ""
    plot(x1, y1, log = logx, xlim = xlim, xlab = xlab, ylab = ylab, 
        main = main, type = "n", ...)
    points(x1, y1, pch = pch1, col = col1, cex = cexp)
    points(x2, y2, pch = pch2, col = col2, cex = cexp)
    temp <- ks.test(x1, x2)
    print(temp)
    if (ifresult) {
        H0 <- "accept"
        if (temp$p.value < 0.05) 
            H0 <- "reject"
        text(locator(1), adj = c(0, 1), paste("Two-Sample Kolmogorov-Smirnov Test", 
            "\nH0: The two data sets are drawn\nfrom the same population\n", 
            "\np-value = ", signif(temp$p.value, 4), "\n", H0, 
            " H0 at 0.05 level", sep = ""), cex = cex)
    }
    label1 <- paste(x1lab, ", N = ", nx1, sep = "")
    label2 <- paste(x2lab, ", N = ", nx2, sep = "")
    legend(locator(1), c(label1, label2), pch = c(pch1, pch2), 
        col = c(col1, col2), bty = "n", cex = cex)
    invisible()
}

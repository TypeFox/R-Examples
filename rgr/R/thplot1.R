thplot1 <-
function (x1, x2, xname = "", ifzero = 0.01, xlow = NA, xhih = NA, 
    yhih = NA, rsd = 5, ptile = 95, main = "", ...) 
{
    if (length(x1) != length(x2)) {
        cat("\nLengths of vectors not the same\n")
        return()
    }
    temp.x <- remove.na(cbind(x1, x2))
    a1 <- temp.x$x[1:temp.x$n, 1]
    a2 <- temp.x$x[1:temp.x$n, 2]
    ndup <- temp.x$n
    xdif <- abs(a1 - a2)
    xdif[xdif <= 0] <- ifzero
    xbar <- (a1 + a2)/2
    if (is.na(xlow)) 
        xlow <- min(xbar)
    if (is.na(xhih)) 
        xhih <- max(xbar)
    if (is.na(yhih)) 
        yhih <- max(xdif)
    if (yhih == ifzero) 
        yhih <- ifzero * 10
    plot(xbar, xdif, xlim = c(xlow, xhih), ylim = c(ifzero, yhih), 
        xlab = "Mean of Duplicates", ylab = "Absolute Difference between Duplicates", 
        log = "xy", type = "n", main = main, ...)
    points(xbar, xdif)
    text(locator(1), paste("No. of", xname, "duplicates =", ndup), 
        cex = 0.8, adj = 0)
    if (rsd <= 0) 
        return()
    calc <- qnorm(1 - (1 - ptile/100)/2) * rsd * 0.014142
    ylcalc <- calc * xlow
    yhcalc <- calc * xhih
    lines(c(xlow, xhih), c(ylcalc, yhcalc))
    ratio <- xbar/xdif
    for (i in 1:ndup) {
        if (ratio[i] <= xlow/ylcalc) 
            ratio[i] <- 1
        else ratio[i] <- 0
    }
    nout <- sum(ratio)
    test <- binom.test(nout, ndup, 1 - (ptile/100), "greater")
    test.prob <- test$p.value
    if(test.prob >= 0.9999) test.prob <- 0.9999
    label <- paste("RSD =", rsd, "% (2SD precision =", 2 * rsd, 
        "%)\nPercentile =", ptile, "%\nDuplicates 'outside' =", 
        nout, "\nProbability =", signif(test.prob, 4))
    text(locator(1), labels = label, cex = 0.8, adj = c(0, 1))
    invisible()
}

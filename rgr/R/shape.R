shape <-
function (xx, xlab = deparse(substitute(xx)), log = FALSE, xlim = NULL, 
    nclass = NULL, ifbw = FALSE, wend = 0.05, ifnright = TRUE, 
    colr = 8, cex = 0.8, ...) 
{
    frame()
    oldpar <- par()
    on.exit(par(oldpar))
    par(mfrow = c(2, 2), pty = "m", cex.main = 0.9)
    temp.x <- remove.na(xx)
    x <- temp.x$x[1:temp.x$n]
    nobs <- temp.x$n
    if ((is.null(nclass)) && (nobs <= 500)) nclass <- "scott"
    if ((is.null(nclass)) && (nobs > 500)) nclass <- "fd"
    save <- gx.hist(x, xlab = xlab, ylab = " ", log = log, xlim = xlim, 
        main = "Histogram", nclass = nclass, ifnright = ifnright, 
        cex = cex, colr = colr, ...)
    xlim <- save$xlim
    if (ifbw) 
        banner <- "Box and Whisker Plot"
    else banner <- "Tukey Boxplot"
    bxplot(x, xlab = xlab, log = log, xlim = xlim, main = banner, 
        ifbw = ifbw, col = colr, wend = wend, cex = cex, colr = colr, 
        ...)
    gx.ecdf(x, xlab = xlab, ylab = " ", log = log, xlim = xlim, 
        main = "Empirical Cumulative Distribution\nFunction (ECDF)", 
        cex = cex, ...)
    cnpplt(x, xlab = xlab, ylab = " ", log = log, xlim = xlim, 
        main = "% Cumulative Percentage\n(Normal) Probability (CPP) Plot", 
        ifshape = TRUE, cex = cex, ...)
    invisible()
}

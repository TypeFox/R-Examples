qqANOVA <- function (x, y, plot.it = TRUE, xlab = deparse(substitute(x)), 
    ylab = deparse(substitute(y)), ...) 
{
    sx <- sort(x)
    sy <- sort(y)
    lenx <- length(sx)
    leny <- length(sy)
    if (leny < lenx) 
        sx <- approx(1L:lenx, sx, n = leny)$y
    if (leny > lenx) 
        sy <- approx(1L:leny, sy, n = lenx)$y
    if (plot.it) 
        palette(rainbow(n=length(y)))
        plot(sx, sy, xlab = xlab, ylab = ylab, 
            col=as.numeric(factor(sy)), cex=3, pch=16, ...)
    invisible(list(x = sx, y = sy))
}

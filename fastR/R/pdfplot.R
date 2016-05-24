#' @export
pdfplot <-
function (pdf = dnorm, xlim = c(0, 1), xlab, ylab, col, lwd = 2, 
    ...) 
{
    if (missing(col)) {
        col = trellis.par.get("add.line")$col
    }
    if (missing(xlab)) {
        xlab = ""
    }
    if (missing(ylab)) {
        ylab = ""
    }
    x <- seq(xlim[1], xlim[2], length = 500)
    xyplot(pdf(x) ~ x, cex = 0.4, xlab = xlab, ylab = ylab, lwd = lwd, 
        ...)
}

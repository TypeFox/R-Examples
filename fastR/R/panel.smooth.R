#' @export
panel.smooth <-
function (x, y, col = trellis.par.get("plot.symbol")$col, col.smooth = trellis.par.get("add.line")$col, 
    bg = NA, pch = trellis.par.get("plot.symbol")$pch, cex = 1, 
    span = 2/3, iter = 3, ...) 
{
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)
}

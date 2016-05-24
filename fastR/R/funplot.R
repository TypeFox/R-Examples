#' @export
funplot <-
function (fun, xlim, ylim, args = list(), lty = 1, lwd = 2, col = trellis.par.get("plot.line")$col, 
    buffer = 0, n = 50, ...) 
{
    if (missing(ylim)) {
        xvals <- seq(xlim[1], xlim[2], length = n)
        yvals = do.call(fun, c(list(x = seq(xlim[1], xlim[2], 
            length = 100)), args))
        ylim <- range(yvals)
        buffer <- buffer * diff(ylim)
        ylim <- ylim + c(-1, 1) * buffer
    }
    xyplot(ylim ~ xlim, panel = function(x, y, ...) {
        panel.mathdensity(fun, args = args, lwd = lwd, lty = lty, 
            col = col, n = n, ...)
    }, ...)
}

plot.noia.linear <-
function (x, loc = 1:x$nloc, effect = TRUE, epistasis = TRUE, 
    ylim = range(GPmap(x)[, 1]) + c(-1, 1) * max(GPmap(x)[, 2]), 
    ...) 
{
    l <- length(loc)
    op <- par(mfrow = c(l, l), mar = c(1.3, 1.5, 1.3, 1))
    for (j in loc) {
        for (i in loc) {
            if (i != j && epistasis) {
                plottwoloci(x, i, j, ylim = ylim, ...)
            }
            if (i == j && effect) {
                plotlocus(x, i, ylim = ylim, ...)
            }
        }
    }
    par(op)
}

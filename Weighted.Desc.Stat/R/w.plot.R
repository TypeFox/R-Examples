w.plot <-
function(x, y, mu, coef.radii, xlim = NULL, ylim = NULL, lwd = NULL, add = NULL, ...)
{
symbols(x, y, mu * coef.radii, inches = FALSE, xlim = xlim, ylim = ylim,lwd = lwd)
}

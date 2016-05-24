eeplot <-
function (x, y, plot.it = TRUE, xlab = deparse(substitute(x)), 
    ylab = deparse(substitute(y)), main = "E-E Plot", ...) 
{
    ey = expectile(y, seq(0.01, 0.99, by = 0.01))
    ex = expectile(x, seq(0.01, 0.99, by = 0.01))
    if (plot.it) {
        plot(ex, ey, main = main, xlab = xlab, ylab = ylab, ...)
        abline(0, 1, col = "red")
    }
    list(x = as.vector(ex), y = as.vector(ey))
}

qddplot <-
function (x, y, remove.ratio = 0.1, differences.range = NA, differences.rangemin = 10, 
    differences.drawzero = TRUE, quantiles.drawhalf = TRUE, quantiles.showaxis = TRUE, 
    line.lwd = 5, xlab = "Quantile", ylab = "Difference", main = "Quantile Differences", 
    ...) 
{
    if (remove.ratio < 0 || remove.ratio >= 0.5) {
        stop("Parameter remove.ratio should be fall into range [0, 0.5).")
    }
    differences.quantile <- calculate.quantiles(x, remove.ratio) - 
        calculate.quantiles(y, remove.ratio)
    if (is.na(differences.range)) {
        differences.range <- calculate.absmax(differences.quantile)
    }
    if (!is.na(differences.rangemin)) {
        if (differences.range < differences.rangemin) {
            differences.range = differences.rangemin
        }
    }
    plot(c(1, 101 - 2 * 100 * remove.ratio), c(-differences.range, 
        differences.range), type = "n", lty = 0, xaxt = "n", 
        xlab = xlab, ylab = ylab, main = main, ...)
    points(differences.quantile, type = "l", lwd = line.lwd, 
        ...)
    if (differences.drawzero) {
        abline(0, 0, lty = "dotted")
    }
    if (quantiles.drawhalf) {
        abline(v = 51 - 100 * remove.ratio, lty = "dashed")
    }
    if (quantiles.showaxis) {
        if (remove.ratio < 0.2) {
            axis(side = 1, at = seq(1 - 100 * remove.ratio, 101 - 
                100 * remove.ratio, 10), labels = seq(0, 1, 0.1))
        }
        else if (remove.ratio < 45) {
            axis(side = 1, at = seq(1 - 100 * remove.ratio, 101 - 
                100 * remove.ratio, 5), labels = seq(0, 1, 0.05))
        }
        else {
            axis(side = 1, at = seq(1 - 100 * remove.ratio, 101 - 
                100 * remove.ratio, 1), labels = seq(0, 1, 0.01))
        }
    }
}

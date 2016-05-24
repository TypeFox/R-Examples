plot.unidimTest <-
function (x, xlab = "Eigenvalue number", ylab = "Eigenvalue", ...) {
    y1 <- x$Tobs
    y2 <- colMeans(x$T.boot, na.rm = TRUE)
    x <- seq_along(y1)
    matplot(x, cbind(y1, y2), xlab = xlab, ylab = ylab, ...)
    invisible()
}

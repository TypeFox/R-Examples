plot.episplineDensity <- function (x, ...) 
{
plot (x$x.pts, x$y.est, type = "l", xlab = "X", ylab = "Density Estimate", main = "Epispline Density Estimate")
points (x$x, rep (0, length (x$x)), col="red", pch=1, cex=1.5)
}

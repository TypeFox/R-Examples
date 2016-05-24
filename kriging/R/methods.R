#
# Ordinary Kriging
#

# New S3 methods for kriging class extends image/plot generic functions
image.kriging <- function(x, main = NULL, xlab = "", ylab = "", col = heat.colors(100), ...) {
  xseq <- round(x[["map"]]$x / x$pixel)
  yseq <- round(x[["map"]]$y / x$pixel)
  a <- min(xseq):max(xseq)
  b <- min(yseq):max(yseq)

  c <- imageloop(xseq, yseq, x[["map"]]$pred, a, b)
  
  image(a*x$pixel, b*x$pixel, c, main=main, xlab=xlab, ylab=ylab, col=col, ...)
}

plot.kriging <- function(x, main = "Semivariogram", xlab = "Distance", ylab = "Semivariance", ...) {
  plot(x[["semivariogram"]]$distance, x[["semivariogram"]]$semivariance, main=main, xlab=xlab, ylab=ylab, xlim=c(0, max(x[["semivariogram"]]$distance)), ylim=c(0, max(x[["semivariogram"]]$semivariance)))

  model <- x[["model"]]
  sill <- x[["sill"]]
  nugget <- x[["nugget"]]
  range <- x[["range"]]
  a <- 1/3

  if(model=="spherical") {
    curve((sill-nugget) * (1.5*(x/range) - 0.5*(x/range)^3) + nugget, add=T)
  }
  if(model=="exponential") {
    curve((sill-nugget) * (1 - exp(-x/(range*a))) + nugget, add=T)
  }
  if(model=="gaussian") {
    curve((sill-nugget) * (1-exp(-(x^2))/((range^2)*a)) * nugget, add=T)
  }
}

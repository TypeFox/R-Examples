# create an image
# some browsers and/or web server are more happy if you call you script foo.png.R instead of foo.R so they can guess the MIME type from the name
run <- function(...) {
  p <- WebPlot(600, 600)
  plot(rnorm(100), rnorm(100), pch = 19, col = 2)
  p
}

plot.wprof <-
function(x,  shape = c("square", "circle", "equispaced"), ...) {
  Z <- getzeta(x)
  plot(Z, shape, ...)
}

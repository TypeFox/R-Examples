plot.incidence <-
function(x, shape = c("square", "circle", "equispaced"), ...) {
  C <- incidence2cover(x)
  plot(C, shape, ...)
}

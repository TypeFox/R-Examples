print.summary.genlasso <- function(x, ...) {
  class(x) = "matrix"
  print(x, digits=4, print.gap=3)
}

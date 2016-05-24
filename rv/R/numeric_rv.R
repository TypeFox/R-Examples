


is.numeric.rv <- function (x) {
  all(rvsimapply(x, is.numeric))
}

as.numeric.rv <- function (x, ...) {
  simapply(x, as.numeric, ...)
}

as.numeric.rvfactor <- function (x, ...) {
  simapply(x, as.numeric, ...)
}


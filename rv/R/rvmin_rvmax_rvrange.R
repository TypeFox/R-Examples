

rvmin <- function (x) {
  rvsimapply(x, min, na.rm=FALSE)
}

rvmax <- function (x) {
  rvsimapply(x, max, na.rm=FALSE)
}

rvrange <- function (x) {
  rvsimapply(x, range, na.rm=TRUE)
}


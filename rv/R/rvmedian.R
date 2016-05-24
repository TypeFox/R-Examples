

rvmedian <- function (x) {
  UseMethod("rvmedian")
}

rvmedian.rv <- function (x) {
  rvsimapply(x, median, na.rm=TRUE)
}

rvmedian.rvsummary <- function (x) {
  rvquantile(x, probs=0.50)
}




mean.rv <- function(x, ...) {
  rvsims(rowMeans(sims(x), ...))
}



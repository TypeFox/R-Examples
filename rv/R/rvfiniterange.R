
rvfiniterange <- function (x) {
  if (is.rvsummary(x)) {
    apply(sims(x), 2, range)
  } else {
    rvsimapply(x, finiterange)
  }
}



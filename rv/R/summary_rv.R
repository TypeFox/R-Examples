
summary.rv <- function (object, ...) {
  summary(as.rvsummary(object), ...)
}

summary.rvfactor <- function (object, all.levels=TRUE, ...) {
  summary(as.rvsummary(object), all.levels=all.levels, ...)
}



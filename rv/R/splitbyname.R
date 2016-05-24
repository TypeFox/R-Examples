
splitbyname <- function (x) {
  a <- split(x, f = .shortnames(x))
  lapply(a, .setDimensionByName)
}


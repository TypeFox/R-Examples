

# this code is a copy the function parallel:::splitRows so that I can pass an R CMD check
# warning.

splitRows <- function(x, ncl) {
  lapply(parallel::splitIndices(nrow(x), ncl), function(i) x[i, , drop = FALSE])
}

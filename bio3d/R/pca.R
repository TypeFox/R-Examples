"pca" <- function(...) {
  dots <- list(...)
  if(inherits(dots[[1]], "matrix")) {
    class(dots[[1]]) <- c("matrix", "xyz")
    UseMethod("pca", dots[[1]])
  }
  UseMethod("pca")
}


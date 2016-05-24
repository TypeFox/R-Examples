# check for packages, and stop if not installed
check4pkg <- function(x) {
  if (!requireNamespace(x, quietly = TRUE)) {
    stop("Please install ", x, call. = FALSE)
  } else {
    invisible(TRUE)
  }
}

ct <- function(l) Filter(Negate(is.null), l)

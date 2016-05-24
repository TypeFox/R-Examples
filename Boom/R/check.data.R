## A collection of utilities for checking whether objects match
## various concepts.

check.scalar.probability <- function(x) {
  okay <- is.numeric(x) && length(x) == 1 && x >= 0 && x <= 1
  if (!okay) {
    stop("expected a scalar probability")
  }
  return(TRUE)
}

check.positive.scalar <- function(x) {
  okay <- is.numeric(x) && length(x) == 1 && x > 0
  if (!okay) {
    stop("expected a positive scalar")
  }
  return(TRUE)
}

check.scalar.integer <- function(x) {
  okay <- is.numeric(x) && (length(x) == 1) && (abs(x - as.integer(x)) < 1e-10)
  if (!okay) {
    stop("")
  }
}

check.nonnegative.scalar <- function(x) {
  okay <- is.numeric(x) && length(x) == 1 && x >= 0
  if (!okay) {
    stop("expected a non-negative scalar")
  }
  return(TRUE)
}

check.probability.distribution <- function(x) {
  okay <- is.numeric(x) && isTRUE(all.equal(sum(x), 1)) && all(x >= 0)
  if (!okay) {
    stop("expected a discrete probability distribution")
  }
  return(TRUE)
}

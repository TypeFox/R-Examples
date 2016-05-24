rMultinom <- function(size, prob) {
  return (rowSums(stats::rmultinom(n = size %/% (2 ^ 31 - 1),
      size = 2 ^ 31 - 1, prob = prob)) + as.vector(stats::rmultinom(n = 1,
      size = size %% (2^31-1), prob = prob)))
}


rvpois <- function (n=1, lambda) {
  rvvapply(stats:::rpois, n.=n, lambda=lambda)
}


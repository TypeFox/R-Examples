

rvbinom <- function (n=1, size, prob) {
  rvvapply(stats:::rbinom, n.=n, size=size, prob=prob)
}


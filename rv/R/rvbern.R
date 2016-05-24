

rvbern <- function (n=1, prob, logical=FALSE) {
  r <- rvvapply(stats:::rbinom, n.=n, size=1, prob=prob)
  if (logical) {
    r <- as.logical(r)
  }
  return(r)
}


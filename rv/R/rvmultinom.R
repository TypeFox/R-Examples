

rvmultinom <- function(n=1, size=1, prob) {
  if (length(prob)<=1) {
    return(rvbinom(n=n, size=size, prob=prob))
  }
  if (anyisrv(n, size, prob)) {
    r <- rvmapply(stats:::rmultinom, n=n, size=size, prob=prob)
  } else {
    n.sims <- getnsims()
    s <- rmultinom(n=n*n.sims, size=size, prob=prob)
    dim(s) <- c(length(s) %/% n.sims, n.sims)
    r <- rvsims(t(s))
    dim(r) <- c(length(prob), n)
    dimnames(r) <- list(names(prob), NULL)
  }
  return(r)
}




rvcat <- function (n=1, prob, levels=NULL) {
  # NAME
  #  rvcat - Sample Categorical Random Variables
  if (anyisrv(n, prob)) {
    x <- rvmultinom(n=n, size=1, prob=prob)
    s <- sims(x, dimensions=TRUE)
    ds <- dim(s)
    s <- as.logical(s)
    dim(s) <- ds
    f <- function (m) row(m)[m]
    a <- apply(s, 1, f)
    r <- if (is.null(dim(a))) rvsims(a) else rvsims(t(a))
  } else {
    n.sims <- getnsims()
    s <- rmultinom(n=n*n.sims, size=1, prob=prob)
    s <- row(s)[as.logical(s)]
    dim(s) <- c(n.sims, length(s) %/% n.sims)
    r <- rvsims(s)
  }
  rvfactor(r, levels=levels)
}

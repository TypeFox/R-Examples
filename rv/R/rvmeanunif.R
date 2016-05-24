

rvmeanunif <- function (n=1, mode=0, scale=1, df) {
  if (n == 1) {
    u <- rvunif(n=df, min=-1, max=1)
    x <- mean(u)
  } else if (n > 1) {
    n.sims <- getnsims()
    M <- array(runif(n=n.sims * n * df, min=-1, max=1), c(n.sims, n, df))
    R <- t(apply(M, MARGIN=1, rowMeans))
    x <- rvsims(R)
  } else {
    stop("n<1")
  }
  return(mode + scale * x)
}

rvtriang <- function (n=1, mode=0, scale=1) {
  rvmeanunif(n=n, mode=mode, scale=scale, df=2)
}


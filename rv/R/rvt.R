

rvt <- function (n = 1, mu = 0, scale = 1, df, ncp, Sigma) { ## CHECK
  if (! missing(Sigma)) {
    t <- .rvmvt(n = n, Sigma = Sigma, df = df)
  } else {
    if (missing(ncp)) {
      t <- rvvapply(stats:::rt, n. = n, df = df)
    } else {
      t <- rvvapply(stats:::rt, n. = n, df = df, ncp=ncp)
    }
    if (scale != 1) 
      t <- (t * scale)
  }
  if (all(mu != 0)) {
    t <- (t + mu) # t + mu, not mu + t (preserves names)
  }
  return(t)
}



.rvmvt <- function (n=1, Sigma, df=1) {
  x <- sqrt(rvchisq(n=n, df=df)/df)
  # DEBUG: But will this work? x is of length n,
  #   but the returned thing is of length n * nrow(Sigma)!
  return(.rvmvnorm(n=n, mean=0, Sigma=Sigma)/x)
}

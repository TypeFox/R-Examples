lnre.goodness.of.fit <- function(model, spc, n.estimated=0, m.max=15)
{
  if (! inherits(model, "lnre")) stop("first argument must belong to a subclass of 'lnre'")
  if (model$multinomial)
    warning("goodness-of-fit test for multinomial sampling not available yet, falling back to Poisson sampling")

  ## automatically reduce number of spectrum elements if variances become too small
  V.Vm <- VVm(model, 1:m.max, N(spc))
  keep <- V.Vm >= 5 # normal approximation to binomial-type distributions should be reasonably good
  if (!all(keep) && missing(m.max)) {
    new.max <- min(which(!keep)) - 1    # first 1:new.max spectrum elements are suitable
    if (new.max < n.estimated + 2) new.max <- n.estimated + 2 # ensure that df >= 3
    m.max <- new.max
  }

  df <- m.max + 1 - n.estimated         # degrees of freedom of the computed chi-squared statistic
  if (df < 1) stop("m.max too small (must be larger than number of estimated parameters)")

  ## compute multivariate chi-squared test statistic, following Baayen (2001), p. 119ff
  N <- N(spc)                           # observed N, V, Vm from spectrum
  V <- V(spc)
  Vm <- Vm(spc, 1:m.max)

  E.Vm <- EVm(model, 1:m.max, N)        # E[Vm(N)]
  E.Vm.2N <- EVm(model, 1:(2*m.max), 2*N)  # E[Vm(2N)]
  E.V <- EV(model, N)                   # E[V(N)]
  E.V.2N <- EV(model, 2 * N)            # E[V(2N)]

  ## construct covariance matrix R
  err.vec <- c(V, Vm) - c(E.V, E.Vm)    # error vector (p. 119)
  cov.Vm.Vk <- diag(E.Vm) -             # "inner part" of covariance matrix (p. 120, (3.61))
    outer(1:m.max, 1:m.max, function (m,k) choose(m+k, m) * E.Vm.2N[m+k] / 2^(m+k))
  cov.Vm.V <- E.Vm.2N[1:m.max] / 2^(1:m.max) # (p. 121, (3.63))
  R <- rbind(c(VV(model, N), cov.Vm.V), # construct covariance matrix R (p. 119, (3.58))
             cbind(cov.Vm.V, cov.Vm.Vk))

  ## compute chi-squared statistic X2 = err.vec %*% R^(-1) %*% err.vec
  x2 <- t(err.vec) %*% solve(R, err.vec)   # (p. 119, (3.59)) -- solve(R, e) returns R^(-1) %*% e, but should be more accurate and stable
  p <- pchisq(x2, df=df, lower.tail=FALSE) # x2 has this chi-squared distribution under H0 (p. 119f, following (3.59))

  result <- data.frame(X2=x2, df=df, p=p)
  rownames(result) <- ""
  result
}

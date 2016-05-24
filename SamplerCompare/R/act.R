# From SamplerCompare, (c) 2010 Madeleine Thompson

# act.R contains functions related to the computation of the
# autocorrelation time of a Markov chain.  See "Graphical Comparison
# of MCMC Samplers" (http://arxiv.org/abs/1011.4457) for a discussion of
# autocorrelation times and examples of how one might use them.

# Computes the autocorrelation time (correlation length) of the
# vector y using an AR(p) model, where p is chosen with AIC.  y
# should be de-meaned before calling this function.  The return value
# is a list with elements:
#
#   act                 autocorrelation time of y
#   act.025, act.975    endpoints of a 95% CI for act
#   se                  standard error for act
#   order               the order of the AR model chosen

ar.act1  <- function(y) {
  # Make sure y is a vector and has at least enough distinct elements
  # to estimate an AR(1) model.

  stopifnot(NCOL(y)==1)
  if (length(unique(y))<5)
    return(list(act=NA, act.025=NA, act.975=NA, se=NA, order=NA))

  # First try arbitrary AIC-chosen order.  If the conditioning on
  # the asymptotic covariance matrix is problematic, choose a smaller
  # maximum order.  The coefficient choice of 0.7 is not principled.

  # NOTE: Is this loop ever repeated?  The condition cutoff used
  # to be much lower, but that was causing ACT underestimates.

  order.max <- NULL
  repeat {
    A <- ar.yw(y, demean=FALSE, order.max=order.max)
    if (A$order==0)   # force an order of at least one
      A <- ar.yw(y, demean=FALSE, order.max=1, aic=FALSE)

    pi <- A$ar
    pi.var <- A$asy.var.coef
    if (kappa(pi.var)<1/sqrt(.Machine$double.eps) || isTRUE(order.max==1))
      break
    order.max <- floor(A$order * 0.7)
  }

  # Use the AR model to compute the autocorrelation time and its
  # CI.  The SE and CI are not equivalent because the CI is asymmetric.

  acf <- matrix(ARMAacf(ar=pi)[2:(A$order+1)])
  act <- (1-sum(pi*acf))/(1-sum(pi))^2

  # Draw from asymptotic distribution of AR coefs.  Uses SVD because
  # pi.var can be ill-conditioned.

  simulation.length <- min(max(40, length(y)), 5000)
  AX <- mvtnorm::rmvnorm(simulation.length, mean=pi, sigma=pi.var, method='svd')
  act.sim <- numeric(simulation.length)

  # Compute simulation.length estimates of the ACT.

  for (i in 1:simulation.length) {
    pi.sim <- AX[i,]
    acf.sim <- ARMAacf(ar=pi.sim)[2:(A$order+1)]
    if (any(abs(polyroot(c(-1,pi.sim)))<1))        # stationary process?
      act.sim[i] <- Inf
    else
      act.sim[i] <- (1-sum(pi.sim*acf.sim))/(1-sum(pi.sim))^2
  }

  # Compute ACT quantiles.

  act.sim[is.na(act.sim)] <- Inf
  act.025 <- as.numeric(quantile(act.sim, 0.025))   # as.numeric to drop name
  act.975 <- as.numeric(quantile(act.sim, 0.975))   # as.numeric to drop name
  se <- (act.975-act.025)/(2*1.96)

  return(list(act=act, se=se, act.025=act.025, act.975=act.975, order=A$order))
}

# See ?ar.act for more information.

ar.act <- function(Y, true.mean=NULL) {

  # Check that parameters are reasonable and estimate mean and
  # variances if they were not specified.

  Y <- as.matrix(Y)
  stopifnot(is.null(true.mean) || ncol(Y)==length(true.mean))

  if (is.null(true.mean))
    mu <- colMeans(Y)
  else
    mu <- true.mean

  # De-mean each column and compute its autocorrelation time with
  # ar.act1.  Return the act, CI, and SE for the slowest-mixing
  # component.

  acts <- sapply(1:ncol(Y), function(i) ar.act1(Y[,i]-mu[i]))
  max.i <- which.max(unlist(acts['act',]))
  if (length(max.i)!=1)   # happens if all ACTs are NA
    return(list(act=NA, act.025=NA, act.975=NA, se=NA, order=NA))
  else
    return(acts[,max.i])
}

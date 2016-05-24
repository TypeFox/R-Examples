# Perform shrinkage on a sample covariance towards a biased covariance
# Author: Brian Lee Yung Rowe
#
# Examples
# 

############################### PUBLIC METHODS ##############################
# Shrink the sample covariance matrix towards the model covariance matrix for
# the given time window.
# model - The covariance matrix specified by the model, e.g. single-index, 
#   Barra, or something else
# sample - The sample covariance matrix. If the sample covariance is null, then
#   it will be computed from the returns matrix
# Example
#   S.hat <- cov_shrink(ys)

cov.shrink(h) %::% TawnyPortfolio : matrix
cov.shrink(h) %as% { cov.shrink(h$returns) }

# This is a general interface for shrinking
cov.shrink(h, prior.fun) %::% a : Function : matrix
cov.shrink(h, prior.fun=cov.prior.cc) %as%
{
  S <- cov.sample(h)

  T <- nrow(h)
  F <- prior.fun(S)
  k <- shrinkage.intensity(h, F, S)
  d <- max(0, min(k/T, 1))

  flog.trace("Got intensity k = %s and coefficient d = %s",k,d)

  S.hat <- d * F + (1 - d) * S
  S.hat
}


# Estimate the covariance matrix using the specified constant.fun for 
# determining the shrinkage constant. The constant.fun is passed two 
# parameters - the sample covariance matrix and the number of rows (T) in the 
# original returns stream.
cov.shrink(m, T, constant.fun, prior.fun) %::% matrix : numeric : Function : Function : matrix
cov.shrink(m, T, constant.fun, prior.fun=cov.prior.cc) %as%
{
  S <- h
  F <- prior.fun(S)
  d <- constant.fun(S, T)

  flog.trace("Got coefficient d = %s",d)

  S.hat <- d * F + (1 - d) * S
  S.hat
}



# Calculate the sample covariance matrix from a returns matrix
# Returns a T x N returns matrix (preferably zoo/xts)
# p.cov <- cov.sample(p)

cov.sample(returns) %::% AssetReturns : matrix
cov.sample(returns) %as%
{
  # X is N x T
  T <- nrow(returns)
  X <- t(returns)
  ones <- rep(1,T)
  S <- (1/T) * X %*% (diag(T) - 1/T * (ones %o% ones) ) %*% t(X)
  Covariance(S)
}

cov.sample(x) %::% a : matrix
cov.sample(x) %as%
{
  cov(x)
}



##############################################################################
# Constant correlation target
# S is sample covariance
cov.prior.cc <- function(S)
{
  r.bar <- cor.mean(S)
  vars <- diag(S) %o% diag(S)
  F <- r.bar * (vars)^0.5
  diag(F) <- diag(S)
  return(F)
}

# This returns a covariance matrix based on the identity (i.e. no correlation)
# S is sample covariance
cov.prior.identity <- function(S)
{
  return(diag(nrow(S)))
}

# Get mean of correlations from covariance matrix
cor.mean <- function(S)
{
  N <- ncol(S)
  cors <- cov2cor(S)
  2 * sum(cors[lower.tri(cors)], na.rm=TRUE) / (N^2 - N)
}

# Calculate the optimal shrinkage intensity constant
# returns : asset returns T x N
# prior : biased estimator
shrinkage.intensity <- function(returns, prior, sample)
{
  p <- shrinkage.p(returns, sample)

  r <- shrinkage.r(returns, sample, p)
  c <- shrinkage.c(prior, sample)
  (p$sum - r) / c
}

# Sum of the asymptotic variances
# returns : T x N (zoo) - Matrix of asset returns
# sample : N x N - Sample covariance matrix
# Used internally.
# S <- cov.sample(ys)
# ys.p <- shrinkage.p(ys, S)
shrinkage.p <- function(returns, sample)
{
  T <- nrow(returns)
  N <- ncol(returns)
  ones <- rep(1,T)
  means <- t(returns) %*% ones / T
  z <- returns - matrix(rep(t(means), T), ncol=N, byrow=TRUE)

  term.1 <- t(z^2) %*% z^2
  term.2 <- 2 * sample * (t(z) %*% z)
  term.3 <- sample^2
  phi.mat <- (term.1 - term.2 + term.3) / T

  phi <- list()
  phi$sum <- sum(phi.mat)
  phi$diags <- diag(phi.mat)
  phi
}



# Estimation for rho when using a constant correlation target
# returns : stock returns
# market : market returns 
# Example
#   S <- cov.sample(ys)
#   ys.p <- shrinkage.p(ys, S)
#   ys.r <- shrinkage.r(ys, S, ys.p)
shrinkage.r <- function(returns, sample, pi.est)
{
  N <- ncol(returns)
  T <- nrow(returns)
  ones <- rep(1,T)
  means <- t(returns) %*% ones / T
  z <- returns - matrix(rep(t(means), T), ncol=N, byrow=TRUE)
  r.bar <- cor.mean(sample)

  # Asymptotic covariance estimator
  term.1 <- t(z^3) %*% z
  term.2 <- diag(sample) * (t(z) %*% z)
  term.3 <- sample * (t(z^2) %*% matrix(rep(1,N*T), ncol=N))
  # This can be simplified to diag(sample) * sample, but this expansion is
  # a bit more explicit in the intent (unless you're an R guru)
  term.4 <- (diag(sample) %o% rep(1,N)) * sample
  script.is <- (term.1 - term.2 - term.3 + term.4) / T

  # Create matrix of quotients
  ratios <- (diag(sample) %o% diag(sample)^-1)^0.5

  # Sum results
  rhos <- 0.5 * r.bar * (ratios * script.is + t(ratios) * t(script.is))

  # Add in sum of diagonals of pi
  sum(pi.est$diags, na.rm=TRUE) 
  + sum(rhos[lower.tri(rhos)], na.rm=TRUE) 
  + sum(rhos[upper.tri(rhos)], na.rm=TRUE)
}

# Misspecification of the model covariance matrix
shrinkage.c <- function(prior, sample)
{
  squares <- (prior - sample)^2
  sum(squares, na.rm=TRUE)
}

##-------------------------------- DEPRECATED -------------------------------##
cov_sample <- function(...) cov.sample(...)

cov_shrink <- function(...) cov.shrink(...)

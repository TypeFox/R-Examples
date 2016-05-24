## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


ar1.llh.C <- function(beta, mu, phi, W, m0, C0, alpha=NULL)
{
  T   = ncol(beta) - 1;
  N.b = length(mu);
  N   = length(m0);
  N.a = N - N.b

  ## Check
  not.ok = rep(0, 8);
  if (not.ok[1] <- length(mu)  != N.b)
    { print("length(mu)!=N.b") ; }
  if (not.ok[2] <- length(phi) != N.b)
    { print("length(phi)!=N.b"); }
  if (not.ok[3] <- (ncol(W) != N.b || nrow(W) != N.b))
    { print("W is not N.b x N.b"); }
  if (not.ok[4] <- length(m0) != N)
    { print("length(m0) != N"); }
  if (not.ok[5] <- (nrow(C0) != N || ncol(C0) != N))
    { print("C0 is not N x N"); }
  if (not.ok[6] <- N.b > N)
    { print("N.b > N"); }
  if (not.ok[7] <- (nrow(beta) != N.b && ncol(beta) != (T+1)))
    { print("beta is not N.b x T"); }

  if (!is.null(alpha)) {
    if (not.ok[8] <- (length(alpha) != N.a))
      { print("length(alpha) != N.a"); }
  }
  else {
    alpha = 0;
  }
  
  if (!prod(!not.ok)) {
    cat("ar1.llh.C: problem.  Returning NA.\n");
    return(NA)
  }    

  log.dens  = 0
  
  OUT <- .C("ar1_llh", alpha, beta,
            mu, phi, W,
            m0, C0,
            as.integer(N.b), as.integer(N), as.integer(T),
            log.dens,
            PACKAGE="BayesLogit");
     
  OUT[[11]]

}

##------------------------------------------------------------------------------

## draw.W.ar1 <- function(x, mu, phi, nu.0, Psi.0)
## {
##   require("bayesm")
##   P = nrow(x)
##   T = ncol(x)
  
##   phi = as.numeric(phi)
##   x = x - as.numeric(mu)
##   e = x[,-1,drop=FALSE] - x[,-T,drop=FALSE] * phi;

##   Psi.T = Psi.0 + e %*% t(e);
##   nu.T  = nu.0 + (T-1);

##   W = rwishart(nu.T, solve(Psi.T))$IW
##   W
## }

draw.W.ar1.ind <- function(x, mu, phi, a.0, b.0)
{
  ## Assume phi is not a matrix.  It needs less structure.
  P = nrow(x)
  T = ncol(x);
  
  phi = as.numeric(phi)
  x = x - as.numeric(mu)
  e = x[,-1,drop=FALSE] - x[,-T,drop=FALSE] * phi;

  a.N = (T-1) + a.0
  b.N = apply(e, 1, function(x){ sum(x^2) }) + b.0;

  W = 1 / rgamma(P, 0.5*a.N, rate=0.5*b.N);
  W
}

draw.phi.ar1 <- function(x, mu, W, phi.m0, phi.P0, phi.prev)
{
  ## Assumes W is square
  P = nrow(x)
  T = ncol(x);
  
  x = x - as.numeric(mu)
  WI = solve(W);
  
  lx.x  = x[,-T,drop=FALSE] %*% t(x[,-1,drop=FALSE]);
  lx.lx = x[,-T,drop=FALSE] %*% t(x[,-T,drop=FALSE]);
  
  b.L = rowSums(lx.x * WI);
  P.L = lx.lx * WI;

  P.N = phi.P0 + P.L
  b.N = phi.P0 %*% phi.m0 + b.L

  U.N = chol(P.N)
  V.N = chol2inv(U.N)
  m.N = V.N %*% b.N;

  phi = m.N + backsolve(U.N, rnorm(P))
  phi

  ok = phi > 0 & phi < 1
  if (!all(ok)) phi = phi.prev

  phi
}

draw.phi.ar1.ind <- function(x, mu, W, phi.m0, phi.P0, phi.prev)
{
  ## Assumes W is diagonal.
  
  ## Dim.
  P = nrow(x)
  N = ncol(x)

  ## To Return.
  phi = 0;

  x = x - as.numeric(mu)
  WI = 1 / W;
  
  lx.x  = apply(x[,-N,drop=FALSE] * x[,-1,drop=FALSE], 1, sum) ## g0
  lx.lx = apply(x[,-N,drop=FALSE] * x[,-N,drop=FALSE], 1, sum) ## g1
  
  b.L = lx.x  * WI;
  P.L = lx.lx * WI;

  P.N = phi.P0 + P.L
  b.N = phi.P0 * phi.m0 + b.L
  V.N = 1 / P.N
  m.N = V.N * b.N
  phi = m.N + sqrt(V.N) * rnorm(P)

  ## Propose and check if accept.
  ok  = phi > 0 & phi < 1
  phi = ifelse (ok, phi, phi.prev);

  phi
}

draw.mu.ar1 <- function(x, phi, W, mu.m0, mu.P0)
{
  ## Assumes W is square.
  P = nrow(x)
  T = ncol(x);

  ## one minus phi
  phi = as.numeric(phi)
  omp = (1-phi)
  WI = solve(W)

  e = x[,-1,drop=FALSE] - x[,-T,drop=FALSE] * phi;
  e = e / omp;
  
  P.L = (T-1) * (omp %*% t(omp)) * WI
  b.L = P.L %*% apply(e, 1, mean);
  
  P.T = P.L + mu.P0
  b.T = b.L + mu.P0 %*% mu.m0;

  U.T = chol(P.T)
  V.T = chol2inv(U.T)
  m.T = V.T %*% b.T;

  mu = m.T + backsolve(U.T, rnorm(P))
  mu
}

draw.mu.ar1.ind <- function(x, phi, W, mu.m0, mu.P0)
{
  ## Assume W is diagonal, i.e. a vector.
  ## Assume phi is not a matrix.
  P = nrow(x)
  T = ncol(x);

  ## one minus phi inverse
  phi = as.numeric(phi)
  omp = (1-phi)
  WI  = 1 / W

  e = x[,-1,drop=FALSE] - x[,-T,drop=FALSE] * phi;
  e = e / omp;
  
  P.L = (T-1) * (omp)^2 * WI
  b.L = P.L * apply(e, 1, mean);

  P.T = P.L + mu.P0
  b.T = b.L + mu.P0 * mu.m0;

  V.T = 1 / P.T
  m.T = V.T * b.T

  mu = m.T + sqrt(V.T) * rnorm(P)
  mu
}

################################################################################
                                   ## TEST ##
################################################################################

if (FALSE) {

  P = 1
  T = 1000

  mu  = rep(2.0, P)
  phi = rep(0.95, P)
  W     = diag(2.0, P)
  W.ind = diag(W)
  W.L = t(chol(W))
  
  m0 = mu
  C0 = diag(1.0, P)

  x = matrix(0, nrow=P, ncol=T+1);
  x[,1] = m0

  for (i in 2:(T+1)) {
    x[,i] = mu + phi * (x[,i-1] - mu) + W.L %*% rnorm(P)
  }

  ## Test:
  samp = 1000

  W.p       = array(0, dim=c(samp, P, P));
  phi.p     = matrix(0, ncol=P, nrow=samp);
  mu.p      = matrix(0, ncol=P, nrow=samp);
  W.p.ind   = matrix(0, ncol=P, nrow=samp);
  phi.p.ind = matrix(0, ncol=P, nrow=samp);
  mu.p.ind  = matrix(0, ncol=P, nrow=samp);

  ## Prior
  mu.m0  = mu
  mu.P0  = diag(10, P)
  phi.m0 = phi
  phi.P0 = diag(10, P)
  W.nu.0 = 0
  W.Psi.0 = diag(0, P)

  mu.m0.ind  = mu
  mu.P0.ind  = rep(10, P)
  phi.m0.ind = phi
  phi.P0.ind = rep(10, P)
  a.0        = rep(1, P)
  b.0        = rep(1, P)

  ## source("AR1.R")
  for (i in 2:samp) {

    mu.p[i,]  = draw.mu.ns(x, phi, W, mu.m0, mu.P0)
    phi.p[i,] = draw.phi.ns(x, mu, W, phi.m0, phi.P0, phi.p[i-1,])
    W.p[i,,]  = draw.W.ns(x, mu, phi, W.nu.0, W.Psi.0)
    
    mu.p.ind[i,]  = draw.mu.ns.ind(x, phi, W.ind, mu.m0.ind, mu.P0.ind)
    phi.p.ind[i,] = draw.phi.ns.ind(x, mu, W.ind, phi.m0.ind, phi.P0.ind, phi.p.ind[i-1,])
    W.p.ind[i,]   = draw.W.ns.ind(x, mu, phi, a.0, b.0)
      
  }

  apply(mu.p, 2, mean)
  apply(phi.p, 2, mean)
  apply(W.p, c(2,3), mean)

  apply(mu.p.ind, 2, mean)
  apply(phi.p.ind, 2, mean)
  apply(W.p.ind, 2, mean)
  
}

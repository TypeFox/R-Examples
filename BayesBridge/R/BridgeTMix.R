######################################################################

## Copyright 2012 Nicholas G. Polson, James G. Scott, Jesse Windle
## Contact info: <jwindle@ices.utexas.edu>.

## This file is part of BayesBridge.

## BayesBridge is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
  
## BayesBridge is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
  
## You should have received a copy of the GNU General Public License
## along with BayesBridge.  If not, see <http:##www.gnu.org/licenses/>.
			      
######################################################################

## Bridge regression using the mixture of triangles formulation.

## library("BayesBridge")
## library("msm");
## library(mvtnorm)
## library(truncnorm)

draw.u <- function(tau, beta, w, alpha)
{
  m = 1-{abs(beta/tau)*w^(-1/alpha)}
  runif(length(beta),max=m)
}

draw.w <- function(tau, beta, u, alpha)
{
  p = length(beta)
  a = (abs(beta/tau)/(1-u))^alpha
  pr = alpha/(1+alpha*a)
  ## p1 = 0.5*(1+alpha)
  ## pr = p1/(1-p1*a)
  shape = (runif(p) < pr) + 1;
  w = rgamma(p, shape, 1);
  w = w+a
  list("w"=w, "shape"=shape);
}

draw.w.marg <- function(alpha, p)
{
  pr = 1-alpha
  shape = (runif(p) < pr) + 1;
  w = rgamma(p, shape, 1)
  list("w"=w, "shape"=shape);
}

draw.tau <- function(beta, alpha, c, d)
{
  p = length(beta)
  nu = rgamma(1, c + p/alpha, rate=d + sum(abs(beta)^alpha))
  tau = nu^(-1/alpha)
  return(tau);
}

draw.sig2 <- function(beta, x, y, sig2.shape=0.0, sig2.scale=0.0)
{
  n = length(y)
  rss = sum( (as.matrix(y)-x%*%as.matrix(beta))^2 )
  prec = rgamma(1, sig2.shape+n/2, rate=sig2.scale+rss/2)
  return(1/prec)
}

##------------------------------------------------------------------------------

sig.for.pg <- function(tau, alpha)
{
  sig2 = tau^2 * gamma(3/alpha) / gamma(1/alpha);
  sqrt(sig2)
}

## mydpgnorm <- function(y, m, tau, alpha, log=FALSE)
## {
##   sig = sig.for.pg(tau, alpha);
##   out = dpgnorm(y, alpha, m, sig);
##   if (log) out = log(out)
##   out
## }

llh.alpha <- function(alpha, s)
{
  p = length(s);
  p * log(alpha) - p * lgamma(1/alpha) - sum(exp(alpha * s))
}

draw.alpha <- function(alpha, beta, tau, pr.a=1.0, pr.b=1.0, ep=0.1)
{
  s = log(abs(beta / tau));
  a.old = alpha;
  lb = 0.0
  ub = 1.0

  ## We might want to change these bounds if we have some belief about the value
  ## of alpha.
  l.new = max(c(lb, a.old-ep))
  r.new = min(c(ub, a.old+ep))
  d.new = r.new - l.new
  a.new = runif(1, l.new, r.new);

  d.old = min(c(ub, a.new+ep)) -  max(c(lb, a.new-ep))

  ## s.new = tau * a.new^(-1/a.new);
  ## s.old = tau * a.old^(-1/a.old);
  ##log.accept = sum(mydpgnorm(beta, 0, tau, a.new, log=TRUE)) - sum(mydpgnorm(beta, 0, tau, a.old, log=TRUE)) +
    log.accept = llh.alpha(a.new, s) - llh.alpha(a.old, s) +
    dbeta(a.new, pr.a, pr.b, log=TRUE) - dbeta(a.old, pr.a, pr.b, log=TRUE) +
    log(d.old) - log(d.new);
  
  if (runif(1) < exp(log.accept)) alpha = a.new

  alpha
}

draw.alpha.2 <- function(alpha, beta, tau, ep=0.1)
{
  s = log(abs(beta / tau));
  a.old = alpha;
  lb = 0.01
  ub = 0.99

  ## We might want to change these bounds if we have some belief about the value
  ## of alpha.
  a.new = runif(1, lb, ub);

  d.old = min(c(ub, a.new+ep)) -  max(c(lb, a.new-ep))

  log.accept = - sum(exp(s*a.new)) + sum(exp(s*a.old))
  
  if (runif(1) < exp(log.accept)) alpha = a.new

  alpha
}

##------------------------------------------------------------------------------
## Sampling beta

## Geweke style Gibbs sampling
draw.beta.1 <- function(beta, bhat, xx, sig2, tau, u, w, alpha)
{
  p = length(bhat)
  b = (1-u)*{w^(1/alpha)}*tau
  for ( i in 1:p )
    {
      m = bhat[i] - crossprod(xx[i,-i],beta[-i]-bhat[-i])/xx[i,i]
      v = sig2/xx[i,i]
      beta[i] = rtnorm(1,m,sqrt(v),-b[i],b[i])
    }
  beta
}

## Rodriguez-Yam style Gibbs sampling - using SVD
draw.beta.2 <- function(beta, a, tV, d, sig2, tau, u, w, alpha)
{
  P = length(beta);

  b = (1-u)*{w^(1/alpha)}*tau
  ## b = (1-u) * tau * exp(log(w) / alpha)
  z = tV %*% beta;

  for (i in 1:P) {
    lmax = -Inf
    rmin =  Inf

    for (j in 1:P) {
      vji = tV[i,j];
      vj  = tV[ ,j];
      ## rji = vj %*% z - vji * z[i];
      rji = vj[-i] %*% z[-i];
      Dif = b[j] - rji
      Sum = b[j] + rji
      left  = ifelse(vji > 0, -Sum, -Dif) / abs(vji);
      right = ifelse(vji > 0,  Dif,  Sum) / abs(vji);
      lmax = max(c(lmax,left))
      rmin = min(c(rmin,right))
    }

    if (d[i]!=0) {
      m = a[i] / (d[i]^2);
      s = sqrt(sig2) / d[i];
      z[i] = rtnorm(1, m, s, lmax, rmin);
    } else {
      cat("d =", d[i], "\n");
      z[i] = runif(1, lmax, rmin);
    }

  }

  beta = t(tV) %*% z;
  beta
}

## Alternate Rodriguez-Yam style Gibbs sampling.  Not as good.
draw.beta.3 <- function(beta, bhat, xx, sig2, tau, u, w, alpha)
{
  P = length(beta)
  evd = eigen(xx)
  rt  = evd$vectors %*% diag(sqrt(evd$values), P) %*% t(evd$vectors);  ## sqrt XX
  irt = evd$vectors %*% diag(sqrt(1/evd$values), P) %*% t(evd$vectors);
  b = (1-u)*{w^(1/alpha)}*tau

  z = rt %*% beta
  m = rt %*% bhat

  for (i in 1:P) {

    left  = rep(0, P)
    right = rep(0, P)

    for (j in 1:P) {
      rji = irt[j,-i] %*% z[-i]
      Dif = b[j] - rji;
      Sum = b[j] + rji;
      left[j]  = ifelse(irt[j,i] > 0, -Sum, -Dif) / abs(irt[j,i]);
      right[j] = ifelse(irt[j,i] > 0,  Dif,  Sum) / abs(irt[j,i]);
    }

    lmax = max(left)
    rmin = min(right)

    z[i] = rtnorm(1, m[i], sqrt(sig2), lmax, rmin);

  }

  beta = irt %*% z;
}

## ## Hamiltonian MCMC
## draw.beta.4 <- function(beta, bhat, xx, sig2, tau, u, w, alpha)
## {
##   require("tmg")
##   p = length(beta);

##   b = (1-u)*{w^(1/alpha)}*tau

##   ## Constraints
##   F = matrix(nrow=2*p, ncol=p)
##   F[1:p,1:p]   = diag(1,p)
##   F[1:p+p,1:p] = -1 * diag(1,p)
##   g = as.vector(c(b,b))

##   ## print(b)
##   ## print(beta)

##   prec = xx / sig2;
##   ell  = prec %*% bhat;

##   out = rtmg(1, prec, ell, beta, F, g, burn.in = 30);

##   drop(out)
## }

################################################################################

bridge.tmix.R <- function(y, X, nsamp, alpha=0.5, sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0,
                          burn=100, sig2=0.0, tau=0.0, verbose=500,
                          beta.true=NULL, omega.true=NULL)
{
  X <- as.matrix(X)
  xx <- t(X)%*%X
  xy <- t(X)%*%y
  ixx <- chol2inv(chol(xx))

  known.sig2 = sig2 > 0
  known.tau  = tau > 0
  known.alpha = alpha > 0
  known.beta = !is.null(beta.true)

  jsvd = svd(X);
  tV = t(jsvd$v);
  d  = jsvd$d
  A  = jsvd$u %*% diag(d);
  a  = t(A) %*% y;

  p <- ncol(X)

  bhat <- drop(ixx%*%xy)
  beta <- rep(1,p)
  beta <- bhat
  ## tau <- 2
  w = rep(2.0,p)
  u = rep(0.5,p)

  if (sig2 <= 0) sig2 = (1/length(y))*sum((y-X%*%bhat)^2)
  if (tau  <= 0) tau  = 1;
  if (alpha <=0) alpha = 0.5;

  if (known.beta) beta = beta.true;

  output <- list(u = matrix(nrow=nsamp, ncol=p),
                 w = matrix(nrow=nsamp, ncol=p),
                 shape = matrix(nrow=nsamp, ncol=p),
                 beta = matrix(nrow=nsamp, ncol=p),
                 sig2 = rep(0, nsamp),
                 tau = rep(0, nsamp),
                 alpha = rep(0, nsamp)
                 )

  start.time = proc.time();

  for( i in 1:(nsamp+burn))
    {
      if( i%%verbose==0 ) cat("iteration ", i, "\n")
      if( i==(burn+1) ) ess.time = proc.time();
      
      if (!known.tau) tau = draw.tau(beta, alpha, nu.shape, nu.rate)
      ## cat("tau", tau, "\n");
      
      if (!known.sig2) sig2 = draw.sig2(beta, X, y, sig2.shape, sig2.scale)
      ## cat("sig2", sig2, "\n");
      
      ws <- draw.w(tau, beta, u, alpha)
      w = ws$w
      shape = ws$shape

      u <- draw.u(tau, beta, w, alpha)

      ## for (k in 1:1) {
      ## u <- draw.u(tau, beta, w, alpha)
      ## beta = draw.beta.1(beta, bhat, xx, sig2, tau, u, w, alpha)
      if (!known.beta) beta = draw.beta.2(beta, a, tV, d, sig2, tau, u, w, alpha)
      ## beta = draw.beta.3(beta, bhat, xx, sig2, tau, u, w, alpha)
      ## beta = draw.beta.4(beta, bhat, xx, sig2, tau, u, w, alpha);
      ## }

      if (!known.alpha) alpha = draw.alpha(alpha, beta, tau);
      ## cat("alpha:", alpha, "\n")
      
      if(i > burn)
        {
          output$u[i-burn,]    = u
          output$w[i-burn,]    = w
          output$shape[i-burn,] = shape
          output$beta[i-burn,] = beta
          output$sig2[i-burn]  = sig2
          output$tau[i-burn]   = tau
          output$alpha[i-burn] = alpha
        }
    }

  end.time = proc.time();
  output$runtime = (end.time - start.time)[1];

  output
}

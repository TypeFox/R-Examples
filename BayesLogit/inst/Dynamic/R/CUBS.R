## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

## CUBS: conjugate updating backward sampling
## Use a MH step to to correct.

binom.solve <- function(rs, fq)
{
  r = rs[1]
  s = rs[2]
  E = digamma(r) - digamma(s)
  V = trigamma(r) + trigamma(s)
  F1 = E - fq[1]
  F2 = V - fq[2]
  out = c("F1"=F1, "F2"=F2)
  out
}

poisson.solve <- function(rs, fq)
{
  F1 = digamma(rs[1]) - log(abs(rs[2])) - fq[1]
  F2 = trigamma(rs[1]) - fq[2]
  out = c("F1"=F1, "F2"=F2)
  out
}

CUBS.R <- function(z, X, n, mu, phi, W, m0, C0, obs=c("binom", "nbinom", "norm"), max.iter=100)
{
  ## When tracking (beta_t, alpha).  It may be the case that there is no alpha.
  ## z_t = x_t (beta_t, alpha_t) + ep_t, ep_t \sim N(0, V_t).
  ## beta_t = mu + phi * (beta_t - mu) + omega_t, omega_t \sim N(0,W).
  ## alpha_t = alpha_{t-1}
  
  ## z : vector of observations (T)
  ## X : design matrix (T x N)
  ## mu : mu (K)
  ## phi : vector (K)
  ## W : covariance MATRIX of innovations of beta (K x K)
  ## m0 : prior mean on (beta_0, alpha_0) (N)
  ## C0 : prior var on (beta_0, alpha_0) (N x N).

  W = as.matrix(W)
  
  T = length(z);
  N.b = ncol(W);
  N   = ncol(X);
  N.a = N - N.b;
  b.idc = 1:N.b + N.a;
  a.idc = 1:N.a;

  with.alpha = N.a > 0;
  
  m = array(m0, dim=c(N, T+1));
  C = array(C0, dim=c(N, N, T+1));
  R = array(0., dim=c(N, N, T+1));
  a = array(0., dim=c(N, T+1));
  rs = array(0., dim=c(2, T));

  beta = array(0, dim=c(N.b, T+1));
  
  d = c( phi, rep(1, N.a) );
  D = diag(d, N);
  mu = c( mu, rep(0, N.a) );
  big.W = matrix(0, N, N); big.W[b.idc, b.idc] = W;
  if (length(n)==1) n = array(n, dim=T);
  
  ## Filter Forward
  for (i in 2:(T+1)) {
    i.l = i-1;

    a[,i]  = d * m[,i-1] + (1-d) * mu;
    R[,,i] = D %*% C[,,i-1] %*% D  + big.W;

    x.i = t(X[i.l,])
    f.i = x.i %*% a[,i];
    q.i = ( x.i %*% R[,,i] %*% X[i.l,] )[1];
    
    rho.i = R[,,i] %*% X[i.l,];
    A.i = rho.i / q.i;

    ## Conjugate update
    if (obs=="binom") { # Binomial
      rs.i = multiroot(binom.solve, start=c(0.1,0.1), maxiter=max.iter, fq=c(f.i, q.i));
      ## cat("f:", f.i, "q:", q.i, "iter:", rs.i$iter, "precis:",
      ##     rs.i$estim.precis, "r:", rs.i$root[1], "s:", rs.i$root[2], "\n");
      rs[,i.l] = rs.i$root
      rstar.i = rs.i$root[1] + z[i.l];
      sstar.i = n[i.l] - z[i.l] + rs.i$root[2];
      fqstar.i = binom.solve(c(rstar.i, sstar.i), c(0, 0));
      ## ## cat("f,q:", f.i, q.i, "root:", rs.i$root, "f.root:", rs.i$f.root, "\n");
    } else if (obs=="norm") { # Gaussian
      qstar.i = (1 / (1 / q.i + 1 / n[i.l]))
      fstar.i = (f.i / q.i + z[i.l] / n[i.l]) * qstar.i
      fqstar.i = c(fstar.i, qstar.i)
      ## cat(c(f.i, q.i), fqstar.i, "\n")
    } else if (obs=="pois") { # Poisson
      r.i = multiroot(function(r,q){trigamma(r)-q}, start=1, q=q.i, maxiter=max.iter);
      s.i = exp(digamma(r.i$root[1]) - f.i)
      rs.i = list(root=c(r.i$root[1], s.i));
      rs[,i.l] = rs.i$root
      fstar.i = digamma(rs.i$root[1] + z[i.l]) - log(rs.i$root[2] + 1)
      qstar.i = trigamma(rs.i$root[1] + z[i.l])
      fqstar.i = c(fstar.i, qstar.i)
    } else if (obs=="nbinom") { # Neg. Binomial
      fhat.i = f.i - log(n[i.l])
      rs.i = multiroot(binom.solve, start=c(0.1,0.1), fq=c(fhat.i, q.i), maxiter=max.iter);
      rs[,i.l] = rs.i$root
      rstar.i = rs.i$root[1] + z[i.l];
      sstar.i = n[i.l] + rs.i$root[2];
      fqstar.i = binom.solve(c(rstar.i, sstar.i), c(0, 0));
      fqstar.i[1] = fqstar.i[1] + log(n[i.l])
    }
    
    m[,i]  = a[,i] + A.i * (fqstar.i[1] - f.i[1]);
    C[,,i] = R[,,i] + rho.i %*% t(rho.i) * ( (fqstar.i[2] / q.i - 1) / q.i );

  }

  ## cat("m[11] (R): \n", m[,11], "\n");
  ## cat("C[11] (R): \n", C[,,11], "\n");
  
  ## Keep track of log density.
  log.dens = 0
  
  ## Backwards sample
  L = t( chol(C[,,T+1]) );
  ## evd = eigen(C[,,T+1]);
  ## Rt = evd$vectors %*% diag(sqrt(evd$values), N) %*% t(evd$vectors);
  ep    = rnorm(N);
  ## ep    = rep(0, N);
  theta = m[,T+1] + L %*% ep;
  ## alpha = ifelse(with.alpha, theta[a.idc], 0);
  if (with.alpha) alpha = theta[a.idc] else alpha = 0
  beta[,T+1] = theta[b.idc];

  log.dens = log.dens - (0.5 * (t(ep) %*% ep) + sum(log(diag(L))));

  for (i in (T+1):2) {

    A.bs = C[b.idc, b.idc, i-1] %*% (solve(R[b.idc, b.idc, i]) * phi);
    V.bs = C[b.idc, b.idc, i-1] - A.bs %*% R[b.idc, b.idc, i] %*% t(A.bs);
    m.bs = m[b.idc, i-1] + A.bs %*% (beta[,i] - a[b.idc,i]);

    print(V.bs)
    
    L  = t(chol(V.bs));
    ep = rnorm(N.b)
    ## ep = rep(0, N.b);
    beta[,i-1] = m.bs + L %*% ep;

    log.dens = log.dens - (0.5 * (t(ep) %*% ep) + sum(log(diag(L))));
    
  }

  ## Need to return log density as well or r,s to calculate log-density.
  out = list("alpha"=alpha, "beta"=beta, "log.dens"=log.dens, "m"=m, "C"=C, "rs"=rs);
  out
}

##------------------------------------------------------------------------------

CUBS.C <- function(z, X, n, mu, phi, W, m0, C0,
                   obs=c("binom", "nbinom", "norm"),
                   eps.rel=1e-8, max.iter=100)
{
  T   = length(z);
  N.b = length(mu);
  N   = ncol(X);
  N.a = N - N.b;

  W = as.matrix(W)

  ## Check
  not.ok = rep(6);
  if (not.ok[1] <- length(mu)  != N.b)
    { cat("length(mu)!=N.b", N.b, "\n") ; }
  if (not.ok[2] <- length(phi) != N.b)
    { cat("length(phi)!=N.b", N.b, "\n"); }
  if (not.ok[3] <- (ncol(W) != N.b || nrow(W) != N.b))
    { cat("W is not N.b x N.b", N.b, "\n"); }
  if (not.ok[4] <- (nrow(X) != T || ncol(X) != N))
    { cat("X is not T x N", T, N, "\n"); }
  if (not.ok[5] <- length(m0) != N)
    { cat("length(m0) != N", N, "\n"); }
  if (not.ok[6] <- (nrow(C0) != N || ncol(C0) != N))
    { cat("C0 is not N x N", N, "\n"); }
  if (not.ok[7] <- length(n) != T)
    { cat("length(n) != T", T, "\n"); }
  if (not.ok[8] <- length(z) != T)
    { cat("length(z) != T", T, "\n"); }
  if (not.ok[9] <- N.b > N)
    { cat("N.b > N", N.b, N, "\n"); }
    
  if (!prod(!not.ok)) {
    cat("CUBS.C: problem.  Returning NA.\n");
    return(NA)
  }    

  alpha = rep(0, max(N.a, 1));
  beta  = array(0, dim=c(N.b, T+1));

  call.name = paste("cubs", obs[1], sep="_");
  log.dens  = 0
  
  OUT <- .C(call.name, alpha, beta,
            z, X, as.double(n), 
            mu, phi, as.double(W),
            m0, C0,
            as.integer(N.b), as.integer(N), as.integer(T),
            log.dens, eps.rel, as.integer(max.iter),
            PACKAGE="BayesLogit");
     
  out = list("alpha"=OUT[[1]], "beta"=OUT[[2]], "log.dens"=OUT[[14]]);
}

################################################################################
                                  ## CHECK ##
################################################################################

if (FALSE) {
  ## 1-D

  T = 200

  phi = 0.95
  W   = 0.1
  mu  = 0
  b.m0  = mu
  b.C0  = matrix(W / (1-phi^2), 1, 1)

  alpha = 1
  X = matrix(1, nrow=T, ncol=2)
  m0 = c(alpha, b.m0)
  C0 = diag(c(2, diag(b.C0)), 2)
  ## m0 = b.m0
  ## C0 = b.C0
  
  beta = rep(alpha, T+1)
  for (i in 1:T)
    beta[i+1] = phi * beta[i] + rnorm(1, 0, sqrt(W))

  V = 0.2 + runif(T, 0, 0.1)
  y = alpha + beta[-1] + rnorm(T, 0, sqrt(V))

  y.bwd = y[T:1]
  X.bwd = X[T:1,,drop=FALSE]

  samp = 1000
  
  beta.fwd = array(0, dim=c(samp, T+1))
  beta.bwd = array(0, dim=c(samp, T+1))
  alpha.fwd = rep(0, samp)
  alpha.bwd = rep(0, samp)

  ## Make sure things look the same going forwards as backwards w/ stationarity.
  for (i in 1:samp) {
    out = CUBS.C(y, X, V, mu, phi, W, m0, C0, "norm");
    beta.fwd[i,] = as.numeric(out$beta)
    alpha.fwd[i] = as.numeric(out$alpha)
    out = CUBS.C(y.bwd, X.bwd, V, mu, phi, W, m0, C0, "norm");
    beta.bwd[i,] = as.numeric(out$beta)
    alpha.bwd[i] = as.numeric(out$alpha)
  }

  fwd.mean = apply(beta.fwd, 2, mean)
  bwd.mean = apply(beta.bwd, 2, mean)

  alpha.fwd.mean = mean(alpha.fwd)
  alpha.bwd.mean = mean(alpha.bwd)

  par(mfrow=c(2,1))
  plot(beta[-1], type="l")
  lines(fwd.mean[-1], col=2)
  lines(rev(bwd.mean[-1]), col=3)
  
  plot(y, type="l")
  lines(fwd.mean[-1] + alpha.fwd.mean, col=2)
  lines(rev(bwd.mean[-1]) + alpha.bwd.mean, col=3)
  
  alpha.fwd.mean
  alpha.bwd.mean
  
}

##------------------------------------------------------------------------------

##------------------------------------------------------------------------------

if (FALSE) {
  ## 2-D

  T = 200
  P = 2

  phi = rep(0.95, P)
  W   = rep(0.1, P)
  mu  = rep(0, P)
  b.m0  = mu
  b.C0  = diag(W / (1-phi^2), P)

  ## phi = rep(1  , P)
  ## W   = rep(0.1, P)
  ## b.m0 = mu
  ## b.C0 = diag(1, P)
  
  xgrid = seq(-1, 1, length.out=T)
  tX = matrix(0, nrow=P, ncol=T);
  freq = c(0, 2, 4, 6)
  ## freq = c(1, 1.1, 1.2, 1.3)
  for (i in 1:P)
    tX[i,] = cos(freq[i] * pi * xgrid);
  ## tX[1,] = rep(c(0,1), T/2)
  ## tX[2,] = rep(c(1,0), T/2)

  ## tX[1,] = rep(1, T)
  ## tX[2,] = rep(0, T)

  X.dyn = t(tX)
  alpha = 0
  ##X = cbind(1, X.dyn)
  ##m0 = c(alpha, b.m0)
  ##C0 = diag(c(2, diag(b.C0)), P+1)
  X = X.dyn
  m0 = b.m0
  C0 = b.C0
  
  beta = matrix(alpha, ncol=T+1, nrow=P)
  for (i in 1:T)
    beta[,i+1] = phi * beta[,i] + rnorm(P, 0, sqrt(W))

  xbeta = beta[,-1] * t(X.dyn)
  
  V = 0.2 + runif(T, 0, 0.1)
  V = V / 1
  psi = alpha + apply(xbeta, 2, sum) 
  y = psi + rnorm(T, 0, sqrt(V))

  y.bwd = y[T:1]
  X.bwd = X[T:1,,drop=FALSE]

  samp = 1000

  beta.fwd = array(0, dim=c(samp, P, T+1))
  beta.bwd = array(0, dim=c(samp, P, T+1))
  beta.3   = array(0, dim=c(samp, 1, T+1))
  alpha.fwd = rep(0, samp)
  alpha.bwd = rep(0, samp)
  alpha.3   = rep(0, samp)

  ## out.1 = FFBS.C(y, X, V, mu, phi, diag(W, P), m0, C0);
  ## out.2 = FFBS.C(y, X[,1,drop=FALSE], V, mu[1], phi[1], diag(W[1], 1), m0[1], diag(C0[1],1));
  
  ## Make sure things look the same going forwards as backwards w/ stationarity.
  for (i in 1:samp) {
    out = CUBS.C(y, X, V, mu, phi, diag(W, P), m0, C0, "norm");
    ## out = FFBS.C(y, X, V, mu, phi, diag(W, P), m0, C0);
    beta.fwd[i,,] = as.numeric(out$beta)
    alpha.fwd[i] = as.numeric(out$alpha)
    out = CUBS.C(y.bwd, X.bwd, V, mu, phi, diag(W,P), m0, C0, "norm");
    ## out = FFBS.C(y.bwd, X.bwd, V, mu, phi, diag(W,P), m0, C0);
    beta.bwd[i,,] = as.numeric(out$beta)
    alpha.bwd[i] = as.numeric(out$alpha)
    out = FFBS.C(y, X[,1,drop=FALSE], V, mu[1], phi[1], diag(W[1], 1), m0[1], diag(C0[1],1));
    ## beta.3[i,,] = out$beta
    ## alpha.3[i]  = out$alpha
  }

  fwd.mean = apply(beta.fwd, c(2,3), mean)
  bwd.mean = apply(beta.bwd, c(2,3), mean)
  ## b3.mean  = apply(beta.3  , c(2,3), mean)
  
  alpha.fwd.mean = mean(alpha.fwd)
  alpha.bwd.mean = mean(alpha.bwd)

  par(mfrow=c(P+1,1))
  for (i in 1:P) {
    plot(beta[i,-1], type="l")
    lines(fwd.mean[i,-1], col=2)
    lines(rev(bwd.mean[i,-1]), col=3)
    ## readline("<ENTER>")
  }

  psi.fwd = apply(fwd.mean[,-1,drop=FALSE] * tX, 2, sum) + alpha.fwd.mean 
  psi.bwd = apply(bwd.mean[,(T:1)+1,drop=FALSE] * tX, 2, sum) + alpha.bwd.mean
  
  plot(y, type="l")
  lines(psi.fwd, col=2)
  lines(psi.bwd, col=3)
  lines(psi, col=4)

  ## alpha.fwd.mean
  ## alpha.bwd.mean
  
}


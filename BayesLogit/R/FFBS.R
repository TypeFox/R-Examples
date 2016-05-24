## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

## FFBS for generalized dynamic linear models.

FFBS.R <- function(z, X, V, mu, phi, W, m0, C0)
{
  ## When tracking (alpha, beta_t).  It may be the case that there is no alpha.
  ## z_t = x_t (alpha_t, beta_t) + ep_t, ep_t \sim N(0, V_t).
  ## beta_t = mu + phi * (beta_t - mu) + omega_t, omega_t \sim N(0,W).
  ## alpha_t = alpha_{t-1}
  
  ## z : vector of observations (T)
  ## X : design matrix (T x N)
  ## mu : mu (K)
  ## phi : vector (K)
  ## W : covariance MATRIX of innovations of beta (K x K)
  ## V: time varying variances (T)
  ## m0 : prior mean on (beta_0, alpha_0) (N)
  ## C0 : prior var on (beta_0, alpha_0) (N x N).

  T = length(z);
  N.b = length(mu);
  N   = ncol(X);
  N.a = N - N.b;
  b.idc = 1:N.b + N.a;
  a.idc = 1:N.a;
  with.alpha = N.a > 0;
  
  m = array(m0, dim=c(N, T+1));
  C = array(C0, dim=c(N, N, T+1));
  R = array(0., dim=c(N, N, T+1));
  a = array(0., dim=c(N, T+1));

  beta = array(0, dim=c(N.b, T+1));
  
  d = c( rep(1, N.a), phi );
  D = diag(d, N);
  mu = c(rep(0, N.a), mu);
  big.W = matrix(0, N, N); big.W[b.idc, b.idc] = W;
  
  ## Filter Forward
  for (i in 2:(T+1)) {
    i.l = i-1;

    a[,i]  = d * m[,i-1] + (1-d) * mu;
    R[,,i] = D %*% C[,,i-1] %*% D  + big.W;

    tF.i = t(X[i.l,])
    f.i  = tF.i %*% a[,i];

    Q    = tF.i %*% R[,,i] %*% t(tF.i) + diag(V[i.l], 1);
    QI   = solve(Q);

    Rtx = R[,,i] %*% X[i.l,];
    ## xRtx = X[i.l,] %*% Rtx;
    ## QI  = diag(nm.p, n.i) - (nm.p / (1/xRtx + sum(nm.p))) %*% t(nm.p);

    e.i = z[i.l] - f.i;

    ## A.i = R[,,i] %*% t(tF.i) %*% QI;
    A.i = Rtx %*% QI;

    ## We could simplify further.
    m[,i] = a[,i] + A.i %*% e.i;
    ## C[,,i] = R[,,i] - A.i %*% Q %*% t(A.i);
    C[,,i] = R[,,i] - Rtx %*% QI %*% t(Rtx);
    
  }

  ## Backward Sample
  L = t( chol(C[,,T+1]) );
  ## evd = eigen(C[,,T+1]);
  ## Rt = evd$vectors %*% diag(sqrt(evd$values), N) %*% t(evd$vectors);
  theta = m[,T+1] + L %*% rnorm(N);
  alpha = ifelse(with.alpha, theta[a.idc], 0);
  beta[,T+1] = theta[b.idc];
  
  for (i in (T+1):2) {

    ## B = C[,,i-1] %*% (solve(R[,,i]) * d);
    ## theta.V = C[,,i-1] - B %*% R[,,i] %*% t(B);
    ## L = t( chol(theta.V[b.idc, b.idc]) );
    ## e = beta[,i] - a[b.idc,i];
    ## beta.m = m[b.idc,i-1] + B[b.idc, b.idc] %*% e;
    ## beta[,i-1] = beta.m + L %*% rnorm(N.b);

    A.bs = C[b.idc, b.idc, i-1] %*% (solve(R[b.idc, b.idc, i]) * phi);
    V.bs = C[b.idc, b.idc, i-1] - A.bs %*% R[b.idc, b.idc, i] %*% t(A.bs);
    m.bs = m[b.idc, i-1] + A.bs %*% (beta[,i] - a[b.idc,i]);

    L = t(chol(V.bs));
    beta[,i-1] = m.bs + L %*% rnorm(N.b);
  }
  
  list("alpha"=alpha, "beta"=beta);
} ## FFBS

##------------------------------------------------------------------------------

FFBS.C <- function(z, X, V, mu, phi, W, m0, C0)
{
  T   = length(z);
  N.b = length(mu);
  N   = ncol(X);
  N.a = N - N.b;

  ## Check
  not.ok = rep(9);
  if (not.ok[1] <- length(mu)  != N.b)
    { cat("length(mu) =", length(mu), "!= N.b =", N.b, "\n") ; }
  if (not.ok[2] <- length(phi) != N.b)
    { cat("length(phi) =", length(phi), "!= N.b =", N.b, "\n"); }
  if (not.ok[3] <- (ncol(W) != N.b || nrow(W) != N.b))
    { cat("dim(W) =", dim(W), "is not N.b x N.b =", N.b, N.b, "\n"); }
  if (not.ok[4] <- (nrow(X) != T || ncol(X) != N))
    { cat("dim(X) =", dim(X), "is not T x N =", T, N, "\n"); }
  if (not.ok[5] <- length(m0) != N)
    { cat("length(m0) =", length(mu), "!= N =", N, "\n"); }
  if (not.ok[6] <- (nrow(C0) != N || ncol(C0) != N))
    { cat("dim(C0) =", dim(C0), "is not N x N =", N, N, "\n"); }
  if (not.ok[7] <- length(V) != T)
    { cat("length(V) =", length(V), "!= T =", T, "\n"); }
  if (not.ok[8] <- length(z) != T)
    { cat("length(z) =", length(z), "!= T =", T, "\n"); }
  if (not.ok[9] <- N.b > N)
    { cat("N.b =", N.b, ">", "N =", N, "\n"); }
    
  if (!prod(!not.ok)) {
    cat("FFBS.C: problem.  Returning NA.\n");
    return(NA)
  }     

  alpha = rep(0, max(N.a, 1));
  beta  = array(0, dim=c(N.b, T+1));
  ldens = 0.0
  
  OUT <- .C("ffbs", alpha, beta,
            z, X, V, 
            mu, phi, W,
            m0, C0,
            as.integer(N.b), as.integer(N), as.integer(T), ldens,
            PACKAGE="BayesLogit");

  ## void ffbs(double *alpha_, double *beta_,
  ##         double *z_, double *X_, double *mu_, double *phi_, double *W_, double *V_,
  ##         double *m0_, double *C0_, int N_b, int N, int T)
     
  out = list("alpha"=OUT[[1]], "beta"=OUT[[2]], "ldens"=OUT[[14]]);
}

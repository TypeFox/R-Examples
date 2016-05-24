## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

################################################################################
                            ## POSTERIOR BY GIBBS ##
################################################################################

## Bayesian logistic regression
##------------------------------------------------------------------------------
logit.R <- function(y, X, n=rep(1, length(y)),
                    m0=rep(0, ncol(X)), P0=matrix(0, nrow=ncol(X), ncol=ncol(X)),
                    samp=1000, burn=500, verbose=500)
{
  ## X: n by p matrix
  ## y: n by 1 vector, avg response
  ## n: n by 1 vector, # of obs at distinct x

  ## Combine data.
  ## new.data = logit.combine(y, X, n);
  ## y = new.data$y;
  ## X = new.data$X;
  ## n = new.data$n;
  ## n.prior = 0.0;

  X = as.matrix(X);
  y = as.numeric(y)

  p = ncol(X)
  N = nrow(X)

  alpha = (y-1/2)*n

  Z = colSums(X*alpha) + P0 %*% m0;
  ## PsiToBeta = solve(t(X) %*% X) %*% t(X);

  w = rep(0,N)
  ## w = w.known;
  beta = rep(0.0, p)

  output <- list(w = matrix(nrow=samp, ncol=N),
                 beta = matrix(nrow=samp, ncol=p)
                 )

  ## c_k = (1:200-1/2)^2 * pi^2 * 4;

  ## Timing
  start.time = proc.time()
  
  ## Sample
  for ( j in 1:(samp+burn) )
  {
    if (j==burn+1) start.ess = proc.time();
    
    ## draw w
    psi = drop(X%*%beta)
    ## Sum of gamma: poor approximation when psi is large!  Causes crash.
    ## w = rpg.gamma(N, n, psi)
    ## Devroye is faster anyway.
    w = rpg.devroye(N, n, psi);

    ## draw beta - Joint Sample.
    PP = t(X) %*% (X * w) + P0;
    ## U = chol(PP);
    ## m = backsolve(U, Z, transpose=TRUE);
    ## m = backsolve(U, m);
    ## beta = m + backsolve(U, rnorm(p))
    S = chol2inv(chol(PP));
    m = S %*% as.vector(Z);
    beta = m + t(chol(S)) %*% rnorm(p);

    # Record if we are past burn-in.
    if (j>burn) {
        output$w[j-burn,] <- w
        output$beta[j-burn,] <- beta
    }

    if (j %% verbose == 0) { print(paste("LogitPG: Iteration", j)); }
  }

  end.time = proc.time()
  output$total.time = end.time - start.time
  output$ess.time   = end.time - start.ess

  ## Add new data to output.
  output$"y" = y;
  output$"X" = X;
  output$"n" = n;

  output
} ## logit.gibbs.R

## Bayesian logistic regression - Normal Prior
##------------------------------------------------------------------------------

## I include this for a fair comparison with LogitFSF.
draw.beta <- function(z, X, w, b.0=NULL, B.0=NULL, P.0=NULL)
{
  ## z: N x 1 outcomes.
  ## X: N x P design matrix.
  ## b.0: prior mean for beta
  ## B.0: prior variance for beta
  ## P.0: prior precision for beta.
  
  ## FS-F use b to denote means and B to denote variances.

  N = nrow(X);
  P = ncol(X);

  if (is.null(b.0)) b.0 = rep(0.0, P);
  if (is.null(P.0)) P.0 = matrix(0.0, P, P);
  if (!is.null(B.0)) P.0 = solve(B.0);

  P.N = t(X) %*% (X * w) + P.0;
  ## S = solve(PC); ## chol2inv works better for larger P?
  S = chol2inv(chol(P.N));
  m = S %*% (as.vector(z) + P.0 %*% b.0);
  beta = m + t(chol(S)) %*% rnorm(P);

} ## draw.beta

logit.gibbs.np.R <- function(y, X, n=rep(1, length(y)),
                             b.0=NULL, B.0 = NULL, P.0 = NULL,
                             samp=1000, burn=500, verbose=500)
{
  ## X: n by p matrix
  ## y: n by 1 vector, avg response
  ## n: n by 1 vector, # of obs at distinct x

  ## DO NOT USE DEFAULT PRIOR
  y.prior=0.5;
  x.prior=colMeans(as.matrix(X));
  n.prior=0.0;
  
  ## Combine data.
  ## new.data = logit.combine(y, X, n, y.prior, x.prior, n.prior);
  ## y = new.data$y;
  ## X = new.data$X;
  ## n = new.data$n;
  ## n.prior = 0.0;

  ## Don't combine.
  X = as.matrix(X);
  y = as.matrix(y);
  
  ## X = as.matrix(X);

  p = ncol(X)
  N = nrow(X)
  
  ## Default prior parameters.
  if (is.null(b.0)) b.0 = rep(0.0, p);
  if (is.null(P.0)) P.0 = matrix(0.0, p, p);
  if (!is.null(B.0)) P.0 = solve(B.0);

  ## Preprocess.
  alpha = drop((y-1/2)*n)
  Z = colSums(X*alpha)

  w = rep(0,N)
  ## w = w.known;
  beta = rep(0.0, p)

  out <- list(w = matrix(nrow=samp, ncol=N),
              beta = matrix(nrow=samp, ncol=p)
              )

  start.time = proc.time()
  
  ## Sample
  for ( j in 1:(samp+burn) )
  {
    if (j==burn+1) start.ess = proc.time();
    ## draw w
    psi = drop(X%*%beta)
    w = rpg.devroye(N, n, psi);

    ## # draw beta - Joint Sample.
    ## PC = t(X) %*% (X * w) + P.0;
    ## ## S = solve(PC); ## chol2inv works better for larger P?
    ## S = chol2inv(chol(PC));
    ## m = S %*% as.vector(Z);
    ## beta = m + t(chol(S)) %*% rnorm(p);
    beta = draw.beta(Z, X, w, b.0=b.0, P.0=P.0);
    
    # Record if we are past burn-in.
    if (j>burn) {
      out$w[j-burn,] <- w
      out$beta[j-burn,] <- beta
    }

    if (j %% verbose == 0) { print(paste("Iteration", j)); }
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess

  ## ## Add new data to output.
  ## output$"y" = y;
  ## output$"X" = X;
  ## output$"n" = n;

  out
} ## logit.gibbs.np.R

################################################################################
                                 ## TESTING ##
################################################################################

#data = read.table("orings.dat",header=TRUE)
#attach(data)
#failure = 2*failure-1
## x = c(53,56,57,63,66,67,68,69, 70,72,73, 75,76,78,79,80,81)
## y = c( 1, 1, 1, 0, 0, 0, 0, 0,3/4, 0, 0,1/2, 0, 0, 0, 0, 0)
## n = c( 1, 1, 1, 1, 1, 3, 1, 1,  4, 1, 1,  2, 2, 1, 1, 1, 1)
## ans = logit.MCMC(100000,cbind(1,x),y,n)

## hist(ans$beta[,1])
## hist(ans$beta[,2])

## mean(ans$beta[,1])
## mean(ans$beta[,2])

################################################################################

## Copyright 2012 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit.

## BayesLogit is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or any later version.

## BayesLogit is distributed in the hope that it will be useful, but WITHOUT ANY
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
## A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

## You should have received a copy of the GNU General Public License along with
## BayesLogit.  If not, see <http://www.gnu.org/licenses/>.

################################################################################

mlogit.R <- function(y, X, n=rep(1,nrow(as.matrix(y))),
                     m.0=array(0, dim=c(ncol(X), ncol(y))),
                     P.0=array(0, dim=c(ncol(X), ncol(X), ncol(y))),
                     samp=1000, burn=500, verbose=1000)
  
{
  ## N - number of trials
  ## J - number of categories
  ## P - number of covariates
  ## y - N x J-1 matrix.  y_{ij} = fraction of outcomes in jth category on trial i.
  ## X - N x P design matrix
  ## n - N x 1 matrix of rolls
  ## Assume beta_J = 0 for identification.

  X = as.matrix(X);
  y = as.matrix(y);
  
  N = nrow(X);
  P = ncol(X);
  J = ncol(y) + 1;

  out = list(
    beta = array(0, dim=c(samp, P, J-1)),
    w    = array(0, dim=c(samp, N, J-1))
    )

  beta = matrix(0, P, J);
  w    = matrix(0, N, J);
  
  ## Precompute.
  kappa = (y - 0.5)*n;

  b.0 = matrix(0, P, J-1);
  for (j in 1:(J-1)) b.0[,j] = P.0[,,j] %*% m.0[,j];

  ## A = rowSums( exp(X %*% beta[,-1]) );
  for (i in 1:(samp+burn)) {
    for (j in 1:(J-1)) {

      ## For now recompute at each iteration.  Try taking out later.
      A = rowSums( exp(X %*% beta[,-j]) );
      
      c.j   = log(A);
      eta.j = X %*% beta[,j] - c.j;

      ## omega.j
      w[,j] = rpg.devroye(N, n, eta.j);
      
     ## beta.j
      PL.j = t(X) %*% (X * w[,j]);
      bL.j = t(X) %*% (kappa[,j] + c.j * w[,j]);

      P1.j = PL.j + P.0[,,j];
      ## Can speed up using Choleksy.
      V1.j = chol2inv(chol(P1.j));
      m1.j = V1.j %*% (bL.j) + b.0[,j];
      
      sqrtV1.j = t(chol(V1.j));
      beta[,j] = m1.j + sqrtV1.j %*% rnorm(P);

      ## A += exp(X %*% beta[,j]) - exp(X %*% beta[,(j+1) %% (J-1)]) 
      
      ## Store
      if (i > burn) {
        out$beta[i-burn,,j] = beta[,j];
        out$w   [i-burn,,j] = w[,j];
      }
    }
    if (i %% verbose == 0) cat("Finished", i, "\n");
  }

  ## Other info
  ## out$psi = array(0, dim=c(1000, N, J-1));
  ## for (i in 1:samp) {
  ##   out$psi[i,,] = X %*% out$beta[i,,];
  ## }
  out$X = X;
  out$y = y;
  out$n = y;

  out
}

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

## Bridge regression using mixture of normals representation.

#####################source("BridgeTMix.R") ## For draw.tau, draw.sig, etc.

bridge.nmix.R <- function(y, X, nsamp, alpha=0.5, sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0,
                          burn=100, sig2=0.0, tau=0.0, verbose=500,
                          beta.true=NULL, lambda.true=NULL)
{
  ## require("copula"); ## No longer required.  Adapted implementation to BayesBridge.

  ## Set up.
  X <- as.matrix(X)
  xx <- t(X)%*%X
  xy <- t(X)%*%y
  ixx <- chol2inv(chol(xx))

  known.sig2 = sig2 > 0
  known.tau  = tau > 0
  known.alpha = alpha > 0
  known.beta = !is.null(beta.true)
  known.lambda = !is.null(lambda.true)

  p = ncol(X)
  n = length(y)

  ## Initialize.
  bhat <- drop(ixx%*%xy)
  beta <- rep(1,p)
  beta = bhat
  lambda = rep(1, p);
  
  if (sig2 <= 0) sig2 = (1/length(y))*sum((y-X%*%bhat)^2)
  if (tau  <= 0) tau  = 1;
  if (alpha<= 0) alpha = 0.5;

  if (known.beta) beta = beta.true
  if (known.lambda) lambda = lambda.true
  
  output <- list(lambda = matrix(nrow=nsamp, ncol=length(beta)),
                 beta   = matrix(nrow=nsamp, ncol=length(beta)),
                 sig2   = rep(0, nsamp),
                 tau    = rep(0, nsamp),
                 alpha  = rep(0, nsamp)
                 )

  colnames(output$beta) = colnames(X);

  start.time = proc.time();
  
  ## GIBBS
  for( i in 1:(nsamp+burn))
    {
      if( i%%verbose==0 ) cat("iteration ", i, "\n")
      if( i==(burn+1) ) ess.time = proc.time();
      
      ## tau -- marginalized draw.
      if (!known.tau) tau = draw.tau(beta, alpha, nu.shape, nu.rate)
      ## cat("tau:", tau, "\n");
      
      ## sig2
      if (!known.sig2) sig2 = draw.sig2(beta, X, y, sig2.shape, sig2.scale)

      ## lambda
      if (!known.lambda) {
        for(j in 1:p)
          lambda[j] = 2 * retstable.ld(0.5 * alpha, 1.0, beta[j]^2 / tau^2);
          ## using copula:
          ## 2 * retstable(0.5 * alpha, 1.0, beta[j]^2 / tau^2, method="LD");
      }

      ## beta
      if (!known.beta) {
        VInv = (xx + diag(lambda * sig2 / tau^2, p));
        V = solve(VInv);
        U = chol(V) * sqrt(sig2);
        m = V %*% xy;
        beta = drop(m + t(U) %*% rnorm(p))
      }
      
      ## alpha
      if (!known.alpha) alpha = draw.alpha(alpha, beta, tau);
      ## cat("alpha:", alpha, "\n");
      
      if(i > burn)
        {
          output$beta[i-burn,]   = beta
          output$lambda[i-burn,] = lambda
          output$sig2[i-burn]    = sig2
          output$tau[i-burn]     = tau
          output$alpha[i-burn]   = alpha
        }
    }

  end.time = proc.time();
  output$runtime = (end.time - start.time)[1];

  output
}

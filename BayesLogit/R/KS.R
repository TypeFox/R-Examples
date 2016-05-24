## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

################################################################################

## Generate random variate from  Kolmogorov-Smirnov distribution.
## Devroye, p. 161 - used as example of alternating series method.
## Alternating series is a clever way to do accept/reject sampling.
## See Devroy, p. 151 for the basic idea.

## This is slightly different than technique he uses for Jacobi paper.

## Here we assume density is written in either case as
## f(x) = c h(x) (1 - a_1(x) + a_2(x) + ...)

## Left case: (index starts at n=0)
## Devroye breaks up sum on p. 161 to get alternating sum.
## a_n = 4x^2/\pi^2 \exp(- (n^2-1) \pi^2 / (8x^2) ),   n odd  (subtract)
## a_n = (n+1)^2 \exp (- ((n+1)^2-1) \pi^2 / (8x^2) ), n even (add)
## ch  = sqrt{2\pi} \pi^2 / (4x^4) \exp(-\pi^2 / (8x^2))

## Right case: (index starts at n=0)
## a_n = (n+1)^2 \exp(-2x^2 ((n+1)^2-1))
## ch  = 8xe^{-2x^2}

## Method for generating proposals come from Lemmas 5.2 and 5.3.

## Just copying Devroye:

draw.left <- function(t.prime)
{
  while (TRUE) {
  
    accept = FALSE
    
    ## Proposal - Truncated Gamma.
    G = 0;
    while (!accept) {
      E = rexp(2);
      E[1] = E[1] / (1 - 0.5 / t.prime);
      E[2] = 2 * E[2];
      G = t.prime + E[1];
      accept = E[1]^2 <= t.prime * E[2] * (G+t.prime);
      if (!accept) accept = G/t.prime - 1 - log(G/t.prime) <= E[2];
    }

    ## Devroye is clever.
    X = pi / sqrt(8 * G);
    W = 0.0
    Z = 0.5 / G;
    P = exp(-G); ## exp(- \pi^2 / (8 X^2) ).
    n = 1;
    Q = 1;
    U = runif(1);
    
    go = TRUE
    while(go) {
      W = W + Z * Q;
      if (U >= W) return(X);
      n = n + 2;
      Q = P^(n^2-1);
      W = W - n^2 * Q;
      go = U >= W;
    }

  }
    
}

draw.right <- function(trunc)
{
  while (TRUE) {
  
    E = rexp(1)
    U = runif(1)
    X = sqrt(trunc^2 + 0.5*E)
    
    W = 0
    n = 1
    Z = exp(-2 * X^2)

    go = TRUE
    while (go) {
      n = n + 1;
      W = W + n^2 * Z^(n^2-1);
      if (U >= W) return(X)
      n = n + 1;
      W = W - n^2 * Z^(n^2-1);
      go = U >= W
    }
    
  }

}

rks.1 <- function()
{
  TRUNC = 0.75;
  t.prime = pi^2 / (8 * TRUNC^2);

  ## 0.373 = pks(0.75) - Differs from PG method.
  p = 0.373;
  q = 1-p;
  
  X = 0;
  if (runif(1) < p/(p+q))
    X = draw.left(t.prime)
  else
    X = draw.right(TRUNC)

  X
}

rks <- function(N=1)
{
  X = rep(0, N)
  for (i in 1:N)
    X[i] = rks.1()
  X
}

pks.1 <- function(x, lower.tail=TRUE, trunc=1000)
{
  i = 1:trunc
  if (x <= 0) return(0);
  a = 2 * (-1)^(i-1) * exp(-2*i^2*x^2);
  p = sum(a);
  if (lower.tail) p = 1 -p;
  p
}

pks <- function(x, lower.tail=TRUE, trunc=1000)
{
  N = length(x);
  for (i in 1:N)
    x[i] = pks.1(x[i], lower.tail, trunc);
  x
}

################################################################################
                                   ## TEST ##
################################################################################

## TEST 1 ##

if (FALSE) {

  ks.samp = rks(100000)

  xgrid = seq(0, 2, 0.01)
  ygrid = pks(xgrid)
  dygrid = diff(ygrid)

  hist(ks.samp, breaks=100, prob=TRUE)
  lines(xgrid[-1], dygrid/0.01, col=2)
  abline(v=0.75, col=3)
  
}

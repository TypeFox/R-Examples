## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

################################################################################
## PG(1.0, Z) - FOLLOWING DEVROYE ##
################################################################################

TRUNC = 0.64
cutoff = 1 / TRUNC;

## pigauss - cumulative distribution function for Inv-Gauss(mu, lambda).
##------------------------------------------------------------------------------
pigauss <- function(x, mu, lambda)
{
  Z = 1.0 / mu;
  b = sqrt(lambda / x) * (x * Z - 1);
  a = -1.0 * sqrt(lambda / x) * (x * Z + 1);
  y = exp(pnorm(b, log.p=TRUE)) + exp(2 * lambda * Z + pnorm(a, log.p=TRUE));
                                        # y2 = 2 * pnorm(-1.0 / sqrt(x));
  y
}

q.and.p <- function(Z)
{
  fz = pi^2 / 8 + Z^2 / 2;
  p = (0.5 * pi) * exp( -1.0 * fz * TRUNC) / fz;
  q = 2 * exp(-1.0 * Z) * pigauss(TRUNC, 1.0/Z, 1.0);

  list("q"=q, "p"=p, "qdivp"=q/p);
}

mass.texpon <- function(Z)
{
  x = TRUNC;
  fz = pi^2 / 8 + Z^2 / 2;
  b = sqrt(1.0 / x) * (x * Z - 1);
  a = -1.0 * sqrt(1.0 / x) * (x * Z + 1);

  x0 = log(fz) + fz * TRUNC;
  xb = x0 - Z + pnorm(b, log.p=TRUE);
  xa = x0 + Z + pnorm(a, log.p=TRUE);

  qdivp = 4 / pi * ( exp(xb) + exp(xa) );

  1.0 / (1.0 + qdivp);
}

mass.detail <- function(Z)
{
  x = TRUNC;
  fz = pi^2 / 8 + Z^2 / 2;
  b = sqrt(1.0 / x) * (x * Z - 1);
  a = -1.0 * sqrt(1.0 / x) * (x * Z + 1);

  x0 = log(fz) + fz * TRUNC;
  xb = x0 - Z + pnorm(b, log.p=TRUE);
  xa = x0 + Z + pnorm(a, log.p=TRUE);

  qdivp = 4 / pi * ( exp(xb) + exp(xa) );

  m = 1.0 / (1.0 + qdivp);
  p = cosh(Z) * 0.5 * pi * exp(-x0);
  q = p * (1/m - 1);
  
  out = list("qdivp"=qdivp, "m"=m, "p"=p, "q"=q, "c"=p+q, "qdivp2"=q/p);

  out
}

## rtigauss - sample from truncated Inv-Gauss(1/abs(Z), 1.0) 1_{(0, TRUNC)}.
##------------------------------------------------------------------------------
rtigauss <- function(Z, R=TRUNC)
{
  Z = abs(Z);
  mu = 1/Z;
  X = R + 1;
  if (mu > R) {
    alpha = 0.0;
    while (runif(1) > alpha) {
      ## X = R + 1
      ## while (X > R) {
      ##     X = 1.0 / rgamma(1, 0.5, rate=0.5);
      ## }
      E = rexp(2)
      while ( E[1]^2 > 2 * E[2] / R) {
        E = rexp(2)
      }
      X = R / (1 + R*E[1])^2
      alpha = exp(-0.5 * Z^2 * X);
    }
  }
  else {
    while (X > R) {
      lambda = 1.0;
      Y = rnorm(1)^2;
      X = mu + 0.5 * mu^2 / lambda * Y -
        0.5 * mu / lambda * sqrt(4 * mu * lambda * Y + (mu * Y)^2);
      if ( runif(1) > mu / (mu + X) ) {
        X = mu^2 / X;
      }
    }
  }
  X;
}

## rigauss - sample from Inv-Gauss(mu, lambda).
##------------------------------------------------------------------------------
rigauss <- function(mu, lambda)
{
  nu = rnorm(1);
  y  = nu^2;
  x  = mu + 0.5 * mu^2 * y / lambda -
    0.5 * mu / lambda * sqrt(4 * mu * lambda * y + (mu*y)^2);
  if (runif(1) > mu / (mu + x)) {
    x = mu^2 / x;
  }
  x
}

## Calculate coefficient n in density of PG(1.0, 0.0), i.e. J* from Devroye.
##------------------------------------------------------------------------------
a.coef <- function(n,x)
{
  if ( x>TRUNC )
    pi * (n+0.5) * exp( -(n+0.5)^2*pi^2*x/2 )
  else
    (2/pi/x)^1.5 * pi * (n+0.5) * exp( -2*(n+0.5)^2/x )
}

## Samples from PG(n=1.0, psi=Z)
## ------------------------------------------------------------------------------
rpg.devroye.1 <- function(Z)
{
  Z = abs(Z) * 0.5;

  ## PG(1,z) = 1/4 J*(1,Z/2)
  fz = pi^2 / 8 + Z^2 / 2;
  ## p = (0.5 * pi) * exp( -1.0 * fz * TRUNC) / fz;
  ## q = 2 * exp(-1.0 * Z) * pigauss(TRUNC, 1.0/Z, 1.0);

  num.trials = 0;
  total.iter = 0;

  while (TRUE)
    {
      num.trials = num.trials + 1;

      if ( runif(1) < mass.texpon(Z) ) {
        ## Truncated Exponential
        X = TRUNC + rexp(1) / fz
      }
      else {
        ## Truncated Inverse Normal
        X = rtigauss(Z)
      }

      ## C = cosh(Z) * exp( -0.5 * Z^2 * X )

      ## Don't need to multiply everything by C, since it cancels in inequality.
      S = a.coef(0,X)
      Y = runif(1)*S
      n = 0

      while (TRUE)
        {
          n = n + 1
          total.iter = total.iter + 1;
          if ( n %% 2 == 1 )
            {
              S = S - a.coef(n,X)
              if ( Y<=S ) break
            }
          else
            {
              S = S + a.coef(n,X)
              if ( Y>S ) break
            }
        }

      if ( Y<=S ) break
    }

  ## 0.25 * X
  list("x"=0.25 * X, "n"=num.trials, "total.iter"=total.iter)
}

## Sample from PG(n, Z) using Devroye-like method.
## n is a natural number and z is a positive real.
##------------------------------------------------------------------------------
rpg.devroye.R <- function(num=1, n=1, z=0.0)
{
  z = array(z, num);
  n = array(n, num);

  total.trials = 0;

  x = rep(0, num);
  for (i in 1:num) {
    x[i] = 0;
    for (j in 1:n[i]) {
      ## x[i] = x[i] + rpg.devroye.1(z[i])
      temp = rpg.devroye.1(z[i]);
      x[i] = x[i] + temp$x;
      total.trials = total.trials + temp$n;
    }
  }
  x
  list("x"=x, "rate"=sum(n)/total.trials)
}

################################################################################
## PG(1.0, Z) - ACCEPT/REJECT ##
################################################################################


rpg.alt.1 <- function(Z)
{
  alpha = 0.0;
  while ( runif(1) > alpha ) {
    X = rpg.devroye.1(0);
    alpha = exp(-0.5 * (Z*0.5)^2 * X);
  }
  X
}

## Sample PG(1.0, Z) using accept/reject.
##------------------------------------------------------------------------------
rpg.alt.R <- function(num=1, Z=0.0)
{
  Z = array(Z, num);
  x = rep(0, num);
  for (i in 1:num) {
    x[i] = rpg.alt.1(Z[i]);
  }
  x
}

################################################################################
                         ## PG(n, Z) - Sum of Gammas ##
################################################################################

## Sample PG(n, z) using sum of Gammas representation.
##------------------------------------------------------------------------------
rpg.gamma.R <-function(num=1, n=1, z=0.0, trunc=200)
{
  w = rep(0, num);
  c.i = (1:200-1/2)^2 * pi^2 * 4
  a.i = c.i + z^2;
  for(i in 1:num){
    w[i] = 2.0 * sum(rgamma(trunc,n)/a.i)
  }
  w
}

################################################################################
                         ## FOR PURPOSES OF TESTING ##
################################################################################

test.igauss <- function(Z, M=100)
{
  x = c(); y = c();
  for (i in 1:M) x[i] = rtigauss(Z);
  for (i in 1:M) {
    draw = TRUNC + 1;
    while (draw > TRUNC) {
      draw = rigauss(1/Z, 1.0);
    }
    y[i] = draw;
  }
  print(paste(mean(x), mean(y), sd(x), sd(y)));

  list("x"=x, "y"=y);
}

if (FALSE) {

  M = 10000;

  x2 = c()
  for ( i in 1:M ) x2[i] = rJstar();
  x2 = x2/4.0

  x3 = c()
  for ( i in 1:M ) x3[i] = rpg(2.0)
  x3 = x3/4.0

  par(mfrow=c(1,2))
  hist(x2, breaks=40, prob=T)
  hist(x3, breaks=40, prob=T)

  print(paste(mean(x2), mean(x3),
              sd(x2), sd(x3)
              ));

}

## test.igauss(20.0, 10000);

################################################################################
## Covergence ##
################################################################################

if (FALSE) {

  M = 100000
  X = matrix(nrow=M, ncol=3)
  for (i in 1:M) {
    X[i,] = as.numeric(rpg.devroye.1(1.37))
  }
  colMeans(X)

}

################################################################################
## Wald's Theorem Stuff ##
################################################################################

int.l.coef <- function(x, n, z)
{
  j = 2 * n + 1;
  mu = j / z;
  lambda = j^2;

  term.1 = 2 * cosh(z) * exp(-1.0 * z * j);
  term.2 = pigauss(x, mu, lambda)

  term.1 * term.2
}

int.r.coef <- function(x, n, z)
{
  j = n + 0.5;

  term.1 = cosh(z) * pi * j;
  term.2 = 0.5 * (z^2 + j^2 * pi^2);
  term.3 = exp(-1 * x * term.2);
  out = term.1 / term.2 * term.3;

  out
}

## x is truncation point.
int.a <- function(x, n, z)
{
  int.l.coef(x,n,z) + int.r.coef(x,n,z);
}

A = int.a(0.64, 0:1000, 0.0)
sum(A)
prob = (A[-1000] - A[-1]) / sum(A);

################################################################################
## Trying to find best split of z for IGauss ##
################################################################################

ig.diff <- function(z) {
  t = 2 / pi
  pigauss(t, 1/z, 1.0) - pigauss(t, 1/z, 1) * exp(-z) / pgamma(1/t, 1/2, rate=1/2, lower.tail=FALSE);
}

## Using rigauss things look okay.

rtig.temp <- function(Z, R=2/pi, N)
{
  Y = 1:N
  num.iter = 0;
  for (i in 1:N){
    X = R + 1;
    alpha = 0.0;
    while (runif(1) > alpha) {
      num.iter = num.iter + 1;
      E = rexp(2)
      while ( E[1]^2 > 2 * E[2] / R) {
        E = rexp(2)
      }
      X = R / (1 + R*E[1])^2
      alpha = exp(-0.5 * Z^2 * X);
    }
    Y[i] = X
  }

  out = list("X"=Y, "a"=N/num.iter)
  out
}

if (FALSE) {
  
  out <- uniroot(ig.diff, c(0,10))
  
}

################################################################################
                     ## Check rejection rate by sampling ##
################################################################################

if (FALSE) {
  Z = 1.37
  outp <- mass.detail(Z);
  samp <- rpg.devroye.R(100000, 1.0, Z);
  
  1/outp$c
  samp$rate
}

################################################################################
                   ## CALCULATING AVERAGE NUMBER OF DRAWS ##
################################################################################

if (FALSE) {

  tk = TRUNC
  ## tk = 2
  
  ## source("PG.R")
  zgrid = seq(0.0, 6, 0.001);
  out   = zgrid
  p = q = zgrid
  len   = length(zgrid)
  ig.accept = zgrid;
  
  for (i in 1:len) {
    Z = zgrid[i]
    fz = pi^2 / 8 + Z^2 / 2;
    p[i] = cosh(Z) * (0.5 * pi) * exp( -1.0 * fz * tk) / fz;
    q[i] = cosh(Z) * 2 * exp(-1.0 * Z) * pigauss(tk, 1.0/Z, 1.0);
    ig.accept[i] = exp(-1.0 * Z) * pigauss(tk, 1.0/Z, 1.0) / pgamma(1/tk, 1/2, 1/2, lower.tail=FALSE)
    out[i] = p[i] + q[i];
  }

  plot(zgrid, out)
  
}

################################################################################
                   ## CALCULATING AVERAGE PROB OF STOPPING ##
################################################################################

if (FALSE) {

  tk = TRUNC
  ## tk = 2
  
  ## source("PG.R")
  n.terms = 16
  lower = 0:n.terms
  upper = 0:n.terms
  int.a = 0:n.terms

  Z = 1.378
  fz = pi^2 / 8 + Z^2 / 2;
  p = cosh(Z) * (0.5 * pi) * exp( -1.0 * fz * tk) / fz;
  q = cosh(Z) * 2 * exp(-1.0 * Z) * pigauss(tk, 1.0/Z, 1.0);

  c.z = p + q;
  
  for (n in 0:n.terms) {
    d.n = 2*n+1
    mu.n = d.n / Z
    lambda.n = d.n^2
    y.n = 0.5 * (Z^2 + (n + 1/2)^2 * pi^2)
    lower[n+1] = 2 * exp(-d.n*Z) * pigauss(tk, mu.n, lambda.n)
    upper[n+1] = 0.5 * pi * d.n / y.n * exp(-y.n * tk)
    int.a[n+1] = lower[n+1] + upper[n+1]
  }

  prob.n = cosh(Z) * (int.a[-n.terms-1] - int.a[-1]) / c.z
  tail   = rev(cumsum(rev(prob.n)))
  
}

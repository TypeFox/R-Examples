## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

## Draw indicators for binomial logistic regression.

##------------------------------------------------------------------------------

draw.indicators.logis.R <- function(z, lambda, nmix)
{
  ## y.u - N x 1 - latent variable y^u in paper.
  ## lambda = X beta

  res = z - log(lambda)
  log.wds = log(nmix$p) - log(nmix$s);

  ## Safer to work on log scale.  Columns correspond to outcome index i!
  log.post = -0.5 * outer(1/nmix$s, res, "*")^2 + log.wds;
  unnrm.post = exp(log.post);

  ## Now sample. 
  r = apply(unnrm.post, 2, function(x){sample.int(n=nmix$N, size=1, prob=x)})
}  ## draw.indicators.logis.R

draw.indicators.logis.C <- function(z, lambda, nmix)
{
  n = length(z);
  r = rep(0, n);
  
  OUT <- .C("draw_indicators_logistic",
            as.integer(r), as.double(z), as.double(lambda), as.integer(n),
            as.double(nmix$p), as.double(nmix$s), as.integer(nmix$N), PACKAGE="BayesLogit")

  OUT[[1]]
} ## draw.indicators.logis.C

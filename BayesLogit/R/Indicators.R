## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

## Draw indicators for negative binomial.

##------------------------------------------------------------------------------

draw.indicators.R <- function(res, nmix)
{
  ## y.u - N x 1 - latent variable y^u in paper.
  ## lambda = X beta

  nc = length(nmix$v)

  nmix$s = sqrt(nmix$v)

  log.wds = log(nmix$p) - log(nmix$s);

  ## Safer to work on log scale.  Columns correspond to outcome index i!  Watch
  ## out for sd in outer product when adapting to other routines.
  log.post = -0.5 * outer(-nmix$m, res, "+")^2 / nmix$v + log.wds;
  unnrm.post = exp(log.post);

  ## Now sample.
  r = apply(unnrm.post, 2, function(x){sample.int(nc, size=1, prob=x)})
}  ## draw.indicators

draw.indicators.C <- function(res, nmix)
{
  n = length(res);
  r = rep(0, n);

  nc = length(nmix$v)

  if (nc != length(nmix$v) || nc != length(nmix$p) || nc==0) {
    cat("draw.indicators.C: problem with dimensions of m,v,p =",
        length(nc), length(nmix$v), length(nmix$p), "\n");
    return(NA)
  }

  ## if (any(!is.finite(res)))
  ##  cat("draw.indicators.C: Residuls have non-finite value.  Dump:\n", res, "\n");

  OUT <- .C("draw_indicators_generic",
            as.integer(r), as.double(res), as.integer(n),
            as.double(nmix$p), as.double(nmix$m), as.double(sqrt(nmix$v)), nc,
            PACKAGE="BayesLogit")

  OUT[[1]]
} ## draw.indicators.C

draw.indicators = draw.indicators.C

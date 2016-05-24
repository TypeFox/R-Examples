## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


draw.omega.mh <- function(omega, beta, llh, y, tX, ntrials, prior.prec, phi, starts, just.max=FALSE, offset=0, type=0)
{
  ## omega: previous innovations values, B x N matrix.
  ## beta:  previous beta values... starting at 1, B x N matrix.
  ## llh:   previous data frame with l0, l1, l2, psi, pl2
  ## y:     response
  ## tX:    transpose of design matrix
  ## ntrials: scalar, number of trials, assuemd constant across response
  ## prior.prec: prior precision of vec(omega), i.e. W values.
  ## phi:   B-vector for AR(1) process.
  ## starts: where the blocks start, indexed using R indexing.
  ## just.max: just move to maximum
  ## offset: N-vector so that psi_i = x_i beta_i + offset_i.

  N = ncol(tX);
  B = nrow(tX);

  ## Go from 1 index notation to 0 index notation.
  starts = starts[!duplicated(starts)]
  starts = starts[order(starts)]
  starts = starts - 1;
  
  if (nrow(omega)!=B || ncol(omega)!=N) {
    cat("Incorrect dimensions for omega.\n");
    return(NA);
  }

  if (nrow(beta)!=B || ncol(beta)!=N) {
    cat("Incorrect dimensions for beta.\n");
    return(NA)
  }

  if (nrow(llh)!=N) {
    cat("Incorrect dimensions for llh.\n");
    return(NA)
  }

  if (starts[1]<0 || starts[length(starts)]>=N) {
    cat("Incorrect starting indices.\n");
    return(NA);
  }

  if (starts[1]!=0) {
    cat("You should probably start the starts at 0.\n");
  }

  ntrials = array(ntrials, dim=N);
  offset  = array(offset, dim=N);

  psi = llh$psi;
  l0  = llh$l0;
  l1  = llh$l1;
  l2  = llh$l2;
  pl2 = llh$pl2;
  psi.dyn = llh$psi.dyn
  psi.stc = llh$psi.stc

  naccept = 0;
  
  OUT <- .C("draw_omega", omega, beta,
            psi.dyn, psi.stc,
            psi, l0, l1, l2, pl2,
            as.double(y), as.double(tX), as.double(ntrials), as.double(offset),
            as.double(prior.prec), as.double(phi),
            as.integer(starts), as.integer(N), as.integer(B), as.integer(length(starts)),
            as.integer(naccept), as.integer(just.max), as.integer(type[1]),
            PACKAGE="BayesLogit");

  out <- list("omega" = OUT[[1]],
              "beta"  = OUT[[2]],
              "llh"   = data.frame("psi.dyn" = OUT[[3]],
                                   "psi.stc" = OUT[[4]],
                                   "psi" = OUT[[5]],
                                   "l0"  = OUT[[6]],
                                   "l1"  = OUT[[7]],
                                   "l2"  = OUT[[8]],
                                   "pl2" = OUT[[9]]),
              "naccept" = OUT[[20]]
              )

  out
}

draw.stc.beta.mh <- function(beta, llh, y, X, ntrials, b0=NULL, P0=NULL, just.max=FALSE, offset=0, type=0)
{
  ## beta:  previous beta values... B-vec.
  ## llh:   previous data frame with l0, l1, l2, psi, pl2
  ## y:     response
  ## tX:    transpose of design matrix
  ## ntrials: scalar, number of trials, assuemd constant across response
  ## m0:    prior mean.
  ## P0:    prior precision.
  ## just.max: just step to maximum.
  ## offset: offset when calculating psi.

  N = nrow(X);
  B = ncol(X);

  if (is.null(b0)) {
    b0 = rep(0, B);
  }

  if (is.null(P0)) {
    P0 = matrix(0, B, B);
  }
  
  if (length(beta)!=B) {
    cat("Incorrect dimensions for beta.\n");
    return(NA)
  }

  if (nrow(llh)!=N) {
    cat("Incorrect dimensions for llh.\n");
    return(NA)
  }

  if (length(b0)!=B) {
    cat("Incorrect dimensions for b0.\n");
    return(NA);
  }

  if (nrow(P0)!=B || ncol(P0)!=B) {
    cat("Incorrect dimensions for P0.\n");
    return(NA);
  }

  ntrials = array(ntrials, dim=N);
  offset = array(offset, dim=N);

  psi = as.double(llh$psi);
  l0  = as.double(llh$l0);
  l1  = as.double(llh$l1);
  l2  = as.double(llh$l2);
  pl2 = as.double(llh$pl2);
  psi.dyn = llh$psi.dyn
  psi.stc = llh$psi.stc

  naccept = 0;
  
  OUT <- .C("draw_stc_beta", beta,
            psi.dyn, psi.stc,
            psi, l0, l1, l2, pl2,
            as.double(y), as.double(X), as.double(ntrials), as.double(offset),
            as.double(b0), as.double(P0),
            as.integer(N), as.integer(B),
            as.integer(naccept), as.integer(just.max), as.integer(type[1]),
            PACKAGE="BayesLogit")

  out <- list("beta"  = OUT[[1]],
              "llh"   = data.frame("psi.dyn" = OUT[[2]],
                                   "psi.stc" = OUT[[3]],
                                   "psi" = OUT[[4]],
                                   "l0"  = OUT[[5]],
                                   "l1"  = OUT[[6]],
                                   "l2"  = OUT[[7]],
                                   "pl2" = OUT[[8]]),
              "naccept" = OUT[[17]]
              )
            
  out
}

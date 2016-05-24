
# PURPOSE: construct beta-prior for rho over -1,1 interval
#-----------------------------------------------------------
# USAGE: out = beta_prior(a1,a2,rvec);
# where:    rvec = grid over rmin,rmax interval, an n x 1 vector
#           a1 = (optional) prior parameter (default = 1.1)
#           a2 = (optional) prior parameter (default = 1.1)
# RETURNS: out = nx1 vector of prior weights for rho
#          in the interval rmin,rmax
#-----------------------------------------------------------
# NOTES: increasing a1,a2 to 1.5,1.5 or 2.0,2.0 increases
#        prior weight placed on zero for rho, and decreases
#        weights on non-zero values for rho
# to see what the prior looks like:
# rvec = -1:0.01:1;
# a1 = 1.1; a2 = 1.1;
# bprior = beta_prior(rvec',a1,a2);
# plot(rvec,bprior);
#-----------------------------------------------------------
# written by:
# James P. LeSage, 4/2003
# University of Toledo
# 2801 W. Bancroft St,
# Toledo, OH 43606
# jlesage@spatial-econometrics.com
#
# Adapted to R by: 
# Miguel Godinho de Matos
# PhD Student at Dept of EPP
# Carnegie Mellon University

beta_prior <- function(rvec,a1=1.01,a2=1.01){
  B          <- beta(a1,a2)
  num        <- (1+rvec)^(a1-1)
  num        <- num * (1-rvec)^(a2-1)
  den        <- 2^(a1+a2-1)
  results    <- (1/B)*num/den
  results[1] <- 0.001
  results[length(results)] <- 0.001
  return( results )
}

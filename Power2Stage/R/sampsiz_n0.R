# ----- helper functions ---------------------------------------------------
# sample size search start values
#
# Sample size for a desired power, large sample approx.
# original 'large sample' sample size formula if diffmthreshold=0
# author D. Labes
# vectorizes properly if diffm is scalar and se a vector or both are vectors
.sampleN00 <- function(alpha=0.05, targetpower=0.8, ltheta1, ltheta2, diffm, 
                       se, steps=2, bk=2, diffmthreshold=0.04)
{
  z1 <- qnorm(1-alpha)
  # within diffmthreshold diffm will be regarded as 0
  # value diffmthreshold=0.04 corresponds roughly to log(0.96)
  # with lower values there are many steps in sample size around between 0.95 and 1.
  # no longer used. Zhangs formula instead
  # vectorized solution
  z2    <- ifelse(abs(diffm)>diffmthreshold, qnorm(targetpower), 
                  qnorm(1-(1-targetpower)/2))
  diffm <- ifelse(abs(diffm)>diffmthreshold, diffm, 0)

  dn <- ifelse(diffm<0, diffm-ltheta1, diffm-ltheta2)
  # return Inf if diffm outside, dn=0 accomplishes that
  dn <- ifelse((diffm-ltheta1)<1.25e-5 | (ltheta2-diffm)<1.25e-5, 0, dn)
  n0 <- (bk/2)*((z1+z2)*(se*sqrt(2)/dn))^2
  
  # round up to next even >=
  # seems Golkowski has used simple round
  n0 <- steps*ceiling(n0/steps)

  return(n0)
  
}

# ----------------------------------------------------------------------------
# 'large sample' sample size with one additional step via t-distribution
# author D. Labes
# bk = design constant, see known.designs()
.sampleN0 <- function(alpha=0.05, targetpower=0.8, ltheta1, ltheta2, diffm, 
                      se, steps=2, bk=2, diffmthreshold=0.04)
{
  # paranoia? swap
  if (ltheta2<ltheta1) {h <- ltheta2; ltheta2 <- ltheta1; ltheta1 <-h}
  
  z1 <- qnorm(1-alpha)
  # value diffmthreshold=0.04 corresponds roughly to log(0.96)
  # with lower values or 0 there are many steps around between 0.95 and 1
  # in sampleN.TOST, no longer used but Zhang's formula instead
  # vectorized solution
  z2    <- ifelse(abs(diffm)>diffmthreshold, qnorm(targetpower), qnorm(1-(1-targetpower)/2))
  diffm <- ifelse(abs(diffm)>diffmthreshold, diffm, 0)
  
  dn <- ifelse(diffm<0, diffm-ltheta1, diffm-ltheta2)
  # return Inf if diffm outside, dn=0 accomplishes that
  dn <- ifelse((diffm-ltheta1)<1.25e-5 | (ltheta2-diffm)<1.25e-5, 0, dn)
  n0 <- (bk/2)*((z1+z2)*(se*sqrt(2)/dn))^2
  
  if (is.finite(n0) & n0>2){
    # make another step with t-distri
    n0 <- ceiling(n0)
    z1 <- qt(1-alpha, df=n0-2)
    z2    <- ifelse(abs(diffm)>diffmthreshold, qt(targetpower, df=n0-2), 
                    qt(1-(1-targetpower)/2, df=n0-2))
    diffm <- ifelse(abs(diffm)>diffmthreshold, diffm, 0)

    dn <- ifelse(diffm<0, diffm-ltheta1, diffm-ltheta2)
    n0 <- (bk/2)*((z1+z2)*(se*sqrt(2)/dn))^2
    
  }
  # make an even multiple of step (=2 in case of 2x2 cross-over)
  n0 <- steps*trunc(n0/steps)
  
  # minimum sample size will be checked outside
  return(n0)
}

# ----------------------------------------------------------------------------
# Paul Zhang (2003)
# A Simple Formula for Sample Size Calculation in Equivalence Studies 
# Journal of Biopharmaceutical Statistics, 13:3, 529-538
.sampleN0_2 <- function(alpha=0.05, targetpower=0.8, ltheta1, ltheta2, diffm, 
                        se, steps=2, bk=2)
{
  # handle unsymmetric limits, Zhang's c0
  c0 <- 0.5*exp(-7.06*(ltheta1+ltheta2)/(ltheta1-ltheta2))
  # Zhang's formula, large sample
  beta <- 1-targetpower
  z1 <- qnorm(1-alpha)
  fz <- ifelse(diffm<0, c0*exp(-7.06*diffm/ltheta1), c0*exp(-7.06*diffm/ltheta2))
  z2 <- abs(qnorm((1-fz)*beta))
  
  dn <- ifelse(diffm<0, diffm-ltheta1, diffm-ltheta2)
  # return Inf if diffm outside, dn=0 accomplishes that
  dn <- ifelse((diffm-ltheta1)<1.25e-5 | (ltheta2-diffm)<1.25e-5, 0, dn)
  n0 <- (bk/2)*((z1+z2)*(se*sqrt(2)/dn))^2
  
  # make an even multiple of step (=2 in case of 2x2 cross-over)
  n0 <- steps*trunc(n0/steps)

  return(n0)
}

# -----------------------------------------------------------------------------
# variant of Zhangs 'large sample' formula with 'smooth' transition from beta/2 
# to beta and a threshold
# vectorizes properly if diffm is scalar and se a vector or both are vectors
.sampleN0_3 <- function(alpha=0.05, targetpower=0.8, ltheta1, ltheta2, diffm, 
                        se, steps=2, bk=2)
{
  # transform to limits symmetric around zero (if they are not)
  locc    <- (ltheta1+ltheta2)/2
  diffm   <- diffm - locc
  ltheta1 <- ltheta1 - locc
  ltheta2 <- -ltheta1
  delta   <- abs((ltheta2-ltheta1)/2)
  
  z1   <- qnorm(1-alpha)
  beta <- 1-targetpower
  
  c  <- abs(diffm/delta)
  # probability for second normal quantil
  # if c<0.2 is in general a good choice needs to be tested
  p2 <- ifelse(c<0.2, 1-(1-0.5*exp(-7.06*c))*beta, 1-beta)
  z2 <- qnorm(p2)
  # difference for denominator
  dn <- ifelse(diffm<0, diffm-ltheta1, diffm-ltheta2)
  # return Inf if diffm outside, dn=0 accomplishes that
  dn <- ifelse((diffm-ltheta1)<1.25e-5 | (ltheta2-diffm)<1.25e-5, 0, dn)
  
  n0 <- (bk/2)*((z1+z2)*(se*sqrt(2)/dn))^2
  # make an even multiple of steps (=2 in case of 2x2 cross-over)
  n0 <- steps*trunc(n0/steps)

  return(n0)
  
}

#----------------------------------------------------------------------------
# Author: dlabes
# 30-Oct-2012, made available as hidden function in V1.2-09 Aug2015
#----------------------------------------------------------------------------
# Calculate the power according to the algo implemented in PASS:
# see Chow S.C., Liu J.P.,
#     "Design and Analysis of Bioavailability
#      and Bioequvalence Studies", Third edition,
#      CRC Press, Chapman & Hall, Boca Raton (2009)
#      Chapter 9.6
# original paper
# Chen, K.W.; Chow, S.C.; and Li, G.
# A Note on Sample Size Determination for Bioequivalence Studies with 
# Higher-Order Crossover Designs.
# J. Pharmacokinetics and Biopharmaceutics, 25 (6), pages 753-765, 1997

# in short: - one degree of freedom less due to carry-over in model
#             for replicate studies
#           - bk=1.1 in case of 2x2x4 also due to carry-over in model
#           - power via shifted central t-approximation
#           - power via 'exact' in case of 2x2
#             but to be consistent we use also shifted as default

# source('C:/Users/dlabes/workspace/PowerTOST/R/DesignHelpers.R')
# source('C:/Users/dlabes/workspace/PowerTOST/R/power.R')
# source('C:/Users/dlabes/workspace/PowerTOST/R/helper.R') 
# source('C:/Users/dlabes/workspace/PowerTOST/R/helper_dp.R')
# source('C:/Users/dlabes/workspace/PowerTOST/R/OwensQ.R')
# source('C:/Users/dlabes/workspace/PowerTOST/R/OwensQOwen.R')

power.PASS <- function(alpha=0.05, logscale=TRUE, theta0, theta1, theta2, 
                       CV, n, design="2x2", method="shifted")
{
  if (missing(CV)) stop("CV must be given!")
  if (missing(n))  stop("Number of subjects n must be given!")
  # check if design is implemented
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design, " unknown!", call.=FALSE)
  
  # get design characteristics
  ades <- .design.props(d.no)
  #degrees of freedom as expression
  dfe  <- .design.df(ades, robust=FALSE)
  #design const.
  bk    <- ades$bk
  if (design=="2x2x4") bk <- 1.1
  # step size = no of sequences
  steps <- ades$steps
  
  #if (design=="2x2" | design=="2x2x2") method="exact"
  # regularize the method giving
  method <- .powerMethod(method)
  
  # handle log-transformation	
  if (logscale) {
    if (missing(theta0)) theta0 <- 0.95
    if (missing(theta1)) theta1 <- 0.8
    if (missing(theta2)) theta2 <- 1/theta1
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    ldiff   <- log(theta0)
    se      <- CV2se(CV)
  } else { # untransformed
    if (missing(theta1)) theta1 <- -0.2
    if (missing(theta0)) theta0 <- 0.05
    if (missing(theta2)) theta2 <- -theta1
    ltheta1 <- theta1
    ltheta2 <- theta2
    ldiff   <- theta0
    se      <- CV
  }
  # handle n = ntotal if scalar else n's of the sequence groups
  if (length(n) == 1) {
    # total n given    
    # for unbalanced designs we divide the ns by ourself
    # to have only small imbalance (function nvec() from Helper_dp.R)
    if(is.finite(n)) n <- nvec(n=n, grps=ades$steps) else n <- rep(Inf, times=steps)
    if (n[1]!=n[length(n)]){
      message("Unbalanced design. n(i)=", paste(n, collapse="/"), " assumed.")
    } 
  } else {
    if (length(n) != ades$steps) {
      stop("Length of n vector must be ", ades$steps, "!")
    }
    if (any(n<1)) stop("All n(i) have to be >0.")
  }
  nc <- sum(1/n)
  n <- sum(n)
  se.fac <- sqrt(ades$bkni * nc)
  
  df <- eval(dfe)
  # in case of replicate designs one degree of freedom lower
  # due to carry over in the model, also for 2x4x2 - Baalam's
  if (design=="2x4x2" | grepl("replicate", ades$name)) {df <- df - 1} 
  if (any(df<1)) stop("n too small. Degrees of freedom <1!")
  
  pow <- .calc.power(alpha, ltheta1, ltheta2, ldiff, sem=se*se.fac, df, 
                     method=method)
  return( pow )
}



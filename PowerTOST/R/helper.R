#----------------------------------------------------------------------------
# Author: dlabes
#----------------------------------------------------------------------------

# helper function to calculate std err from CV of lognormal data
CV2se <- function(CV) sqrt(log(1.0 + CV^2))
# reverse to CV2se helper function
se2CV <- function(se) sqrt(exp(se*se)-1)
# helper function to calculate mse from CV of lognormal data
CV2mse <- function(CV) log(1.0 + CV^2)
mse2CV <- function(mse) sqrt(exp(mse)-1)

#----------------------------------------------------------------------------
# function to calculate confidence limits of given CV 
#----------------------------------------------------------------------------
CVCL <- function(CV, df, side=c("upper", "lower","2-sided"), alpha=0.05)
{
  ssintra <- log(1.0 + CV^2)*df  #s2*df
  side    <- match.arg(side)
  
  limits <- switch(EXPR=side,
      upper= c(0, ssintra/qchisq(alpha,df)),
      lower= c(ssintra/qchisq(1-alpha,df), Inf),
      c(ssintra/qchisq(1-alpha/2,df), ssintra/qchisq(alpha/2,df)))
  limits <-(sqrt(exp(limits)-1))
  names(limits) <- c("lower CL", "upper CL")
  return(limits)
}

#----------------------------------------------------------------------
# Helper function to calculate CV(T) and CV(R) from a 'pooled' CV(T/R)
# assuming a ratio of the intra-subject variances
# Author: dlabes
#----------------------------------------------------------------------

CVp2CV <- function(CV, ratio=1.5)
{
  if(any(ratio<=0)) stop("ratio(s) must be >0!")
  s2p <- CV2mse(CV) 
  # s2d = (s2WT + s2WR)/2 with s2WT/s2wR=ratio
  # s2d = (ratio*s2WR + s2WR)/2
  # s2d = (ratio+1)*s2WR/2
  s2WR <- s2p*2.0/(ratio+1)
  s2WT <- ratio*s2WR
  # return the vector of CVs
  r <- matrix(mse2CV(c(s2WT, s2WR)), ncol=2)
  colnames(r) <- c("CVwT", "CVwR")
  #rownames(r) <- ratio
  if (nrow(r)<2) r <- as.vector(r)
  r
}

#----------------------------------------------------------------------
# Densitiy of inverse gamma distriubtion
# Adapted dinvgamma() from R package MCMCpack (1.3-3)
#----------------------------------------------------------------------
my_dinvgamma <- function(x, shape, scale = 1) {
  if(shape <= 0 | scale <= 0)
    stop("Shape or scale parameter negative in my_dinvgamma().\n")
  exp(shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - (scale/x))
}

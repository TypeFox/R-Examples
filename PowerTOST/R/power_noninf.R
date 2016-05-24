#-----------------------------------------------------------------------------
# Power calculations based on non-inferiority t-test
# 
# Author: dlabes
#-----------------------------------------------------------------------------

# --------------------------------------------------------------------------
# internal functions:
# power function (working horse)
.power.noninf <- function(alpha, lmargin, diffm, sem, df)
{
  tval <- qt(1-alpha, df)
  # the original abs() function has the effect that in case of diffm<lmargin
  # if lmargin<0 the power of inferiority! is calculated
  tau  <- (diffm-lmargin)/sem
  # in case of diffm=lmargin and se=0
  # tau has the value NaN
  tau[is.nan(tau)] <- 0
  
  if (lmargin>0) tau <- -tau
  return(1 - pt(tval, df, tau))
}

# --------------------------------------------------------------------------
# Power function for non-inferiority t-test (OOST)
power.noninf <- function(alpha=0.025,  logscale=TRUE, margin, theta0, CV, n, 
                         design="2x2", robust=FALSE)
{
  # check if design is implemented
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design, " unknown!", call.=FALSE)
  
  # design characteristics
  ades <- .design.props(d.no)
  #degrees of freedom as expression
  dfe  <- .design.df(ades, robust=robust)
  # design constant
  #bk <- ades$bk    # we use always bkni

  # we use always bkni
  #bk   <- ades$bk
  if (length(n) == 1) {
    # total n given    
    # for unbalanced designs we divide the ns by ourself
    # to have only small imbalance (function nvec() from Helper_dp.R)
    n <- nvec(n=n, grps=ades$steps)
    if (n[1]!=n[length(n)]){
      message("Unbalanced design. n(i)=", paste(n, collapse="/"), " assumed.")
    } 
  } else {
    if (length(n) != ades$steps) {
      stop("Length of n vector must be ", ades$steps, "!")
    }
  }
  
  nc <- sum(1/n)
  n <- sum(n)
  se.fac <- sqrt(ades$bkni * nc)
  
  df   <- eval(dfe)    
  if (any(df<1)) stop("n too small. Degrees of freedom <1!")
  
  # handle log-transformation
  if (logscale) {
    if (missing(margin)) margin <- 0.8
    if (missing(theta0)) theta0 <- 0.95
    # further check of input
    if(margin<=0) stop("With logscale=TRUE margin must be ratio >0")
    if(theta0<=0) stop("With logscale=TRUE theta0 must be ratio >0")
    lmargin <- log(margin)
    diffm   <- log(theta0)
    sedm    <- CV2se(CV)*se.fac
  } else {
    if (missing(margin)) margin <- -0.2
    if (missing(theta0)) theta0 <- -0.05
    lmargin <- margin
    diffm   <- theta0
    sedm    <- CV*se.fac
  }
  return(.power.noninf(alpha=alpha, lmargin=lmargin, diffm=diffm, sem=sedm, df=df))
}

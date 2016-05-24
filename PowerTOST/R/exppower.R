#------------------------------------------------------------------------------
# Authors: dlabes, Benjamin Lang
#------------------------------------------------------------------------------

# --- Approximate "expected" power according to Julious book
# taking into account the uncertainty of an estimated se with 
# dfse degrees of freedom
# Only for log-transformed data.
# Raw function: see the call in .exppower.TOST()
.approx.exppower.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, sem, 
                                  dfse, df) 
{
  tval <- qt(1 - alpha, df, lower.tail = TRUE)
  d1   <- sqrt((diffm-ltheta1)^2/sem^2)
  d2   <- sqrt((diffm-ltheta2)^2/sem^2)
  # in case of diffm=ltheta1 or =ltheta2 and se=0
  # d1 or d2 have then value NaN (0/0)
  d1[is.nan(d1)] <- 0
  d2[is.nan(d2)] <- 0
  pow <- pt(d1,dfse,tval) + pt(d2,dfse,tval) - 1
  # approx. may led to negative expected power
  pow[pow<0] <- 0
  pow
}

# --- Exact implementation of expected power
# Author B.Lang
.exact.exppower.TOST <- function(alpha=0.05, ltheta1, ltheta2, ldiff, se.fac,
                                 se, dfCV, df) 
{
  # infinite df of CV should give expected power identical to (conditional) power   
  if (!is.finite(dfCV)) {
    return(.power.TOST(alpha, ltheta1, ltheta2, ldiff, se.fac*se, df))
  }
  # Define assurance function (expected power)
  f <- function(v) {
    .power.TOST(alpha, ltheta1, ltheta2, ldiff, se.fac*sqrt(v), df) * 
      my_dinvgamma(x = v, shape = dfCV/2, scale = dfCV/2 * se^2)
  }
  len <- if (dfCV < 1000) 0.2 else if (dfCV < 1e+04) 0.05 else 0.01
  l <- se^2 - len
  pow <- ifelse(l <= 0, 0, integrate(f, 0, l)$value) + 
    integrate(f, max(l, 0), l + 2*len)$value + 
    integrate(f, l + 2*len, Inf)$value
  
  pow
}

# --- working horse for exppower.TOST() and expsampleN.TOST()
.exppower.TOST <- function(alpha=0.05, ltheta1, ltheta2, ldiff, se.fac,
                           se, dfCV, df, method="exact")
{
  if (method == "exact") {
    return(.exact.exppower.TOST(alpha, ltheta1, ltheta2, ldiff, se.fac, se,
                                dfCV, df))
  } else if (method == "approx") {
    return(.approx.exppower.TOST(alpha, ltheta1, ltheta2, ldiff, sem=se.fac*se,
                                 dfse=dfCV, df=df))
  } else {
    # this is paranoia since method should be checked by the high level function
    stop("Method '", method, "' unknown!\n", call. = TRUE)
  }
}

# ----------------------------------------------------------------------------
# Main function for expected power
# ----------------------------------------------------------------------------
exppower.TOST <- function(alpha=0.05, logscale=TRUE, theta0, theta1, theta2, 
                          CV, dfCV, n, design="2x2", robust=FALSE, 
                          method=c("exact", "approx")) 
{
  # Check if design is implemented
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design, " unknown!", call.=FALSE)
  
  # Design characteristics
  ades <- .design.props(d.no)
  #df as expression
  dfe  <- .design.df(ades, robust=robust)
  #design const.
  #bk   <- ades$bk # we use always bkni
  if (missing(CV) | missing(dfCV)) 
    stop("CV and df must be given!", call.=FALSE)
  if (missing(n)) 
    stop("Number of subjects must be given!", call.=FALSE)
  if (logscale) {
    if (missing(theta0)) theta0 <- 0.95
    if (missing(theta1)) theta1 <- 0.8 
    if (missing(theta2)) theta2 <- 1/theta1
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    ldiff   <- log(theta0)
    se      <- CV2se(CV)
  } else {
    if (missing(theta0)) theta0 <- 0.05
    if (missing(theta1)) theta1 <- -0.2 
    if (missing(theta2)) theta2 <- -theta1
    ltheta1 <- theta1
    ltheta2 <- theta2
    ldiff   <- theta0
    se      <- CV
#    message("It is assumed that CV is the standard deviation.")
  }
  if (length(n) == 1) {
    # total n given    
    # for unbalanced designs we divide the ns by ourself
    # to have only small imbalance (function nvec() from Helper_dp.R)
    n <- nvec(n=n, grps=ades$steps)
    if (n[1] != n[length(n)]) {
      message("Unbalanced design. n(i)=", paste(n, collapse="/"), " assumed.")
    }
  } else {
    if (length(n) != ades$steps) {
      stop("Length of n vector must be ", ades$steps, "!")
    }
  }
  nc <- sum(1/n)
  n  <- sum(n)
  se.fac <- sqrt(ades$bkni * nc)
  df <- eval(dfe)
  if (any(df < 1)) 
    stop("n too low. Degrees of freedom <1!", call.=FALSE)
  
  # check method
  method <- tolower(method)
  method <- match.arg(method)
  
  pwr <- .exppower.TOST(alpha, ltheta1, ltheta2, ldiff, se.fac, se, dfCV, 
                        df, method)
  pwr # return power
}
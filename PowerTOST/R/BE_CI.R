# source("./R/DesignHelpers.R")
# ------------------------------------------------------------------
# 1-2*alpha confidence interval given point est., CV and n
# ------------------------------------------------------------------
CI.BE <- function(alpha=0.05, pe, CV, n, design="2x2", robust=FALSE)
{
  # does it vectorize propperly?
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design, " unknown!", call.=FALSE)
  # design characteristics
  ades <- .design.props(d.no)
  #df as expression
  dfe  <- .design.df(ades, robust=robust)
  
  mse  <- CV2mse(CV)
  # check if n is ntotal
  if (length(n)==1) {
    df <- eval(dfe)
    hw <- qt(1-alpha, df)*sqrt(mse*ades$bk/n)
  } else {
    # or n per sequence group 
    # check the correct length due to design
    # must be length==number of (sequence) groups, that is coded in steps
    if (length(n)!= ades$steps){
      stop("Length of n vector must be ",ades$steps,"!")
    }
    nc <- sum(1/n)
    n <- sum(n)
    df <- eval(dfe)
    hw <- qt(1-alpha, df)*sqrt(mse*ades$bkni*nc)
  }
  upper <- log(pe) + hw
  lower <- log(pe) - hw
  CI    <- cbind(lower,upper)    # gives a matrix
  CI    <- exp(CI)
  if (nrow(CI)==1) CI <- CI[1,] # gives a vector
  return(CI)
}


# ----------------------------------------------------------------------------
# Fieller CI for values in the original scale
# ----------------------------------------------------------------------------
CI.RatioF <- function(alpha=0.025, pe, CV, CVb, n, design=c("2x2","parallel"))
{
  # check design
  design <- match.arg(design)
  
  if (missing(pe)) stop("point est. (ratio T/R) must be given!")
  if (missing(CV)) stop("CV must be given!")
  if (tolower(design)!="parallel"){
    if (missing(CVb)) stop("CVb(etween) must be given in case of cross-over!",
                           call.=FALSE)
  }
  if (missing(n))  stop("Number of subjects must be given!", call.=FALSE)
  if (n<=2) stop("Number of subjects must be >2!", call.=FALSE)
  
  if (length(n)<2){
    #TODO: check even?
    n1 <- trunc(n/2)
    n2 <- n-n1
  } else {
    n1 <- n[1]
    n2 <- n[1]
  }
  df <- n1+n2-2
  tq <- qt(1-alpha, df)
  if (design=="parallel"){
    # Hauschke formulas with YR=1, then YT=pe
    aT <- CV^2*tq^2/n1
    aR <- CV^2*tq^2/n2
    lower <- (pe-sqrt(aR*pe^2+aT-aT*aR))/(1-aR)
    upper <- (pe+sqrt(aR*pe^2+aT-aT*aR))/(1-aR)
  } else {
    # 2x2 crossover
    # again Hauschke formulas with YR=1, then YT=pe
    s2T <- s2R <- CV^2+CVb^2
    aT  <- aR <- 0.25*(1/n1+1/n2)*s2T*tq^2
    aTR <- 0.25*(1/n1+1/n2)*tq^2*CVb^2
    lower <- (pe-aTR-sqrt((pe-aTR)^2-(pe^2-aT)*(1-aR)))/(1-aR)
    upper <- (pe-aTR+sqrt((pe-aTR)^2-(pe^2-aT)*(1-aR)))/(1-aR)
  }
  # if negative values are obtained the condition meanR is >0 
  # is statistically not true. then no Fieller CI is obtainable
  lower <- ifelse(lower<0, NA, lower)
  upper <- ifelse(upper<0, NA, upper)
  if (any(is.na(lower)) | any(is.na(upper))){
    warning("Confidence set(s) unbounded.")
  }
  CI <- cbind(lower,upper)
  if (nrow(CI)==1) CI <- CI[1,] # gives a vector
  return(CI)
}
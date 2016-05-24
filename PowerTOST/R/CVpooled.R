#------------------------------------------------------------------------------
# Calculations for the pooled CV, weighted by df
# 
# Author: dlabes Jan 2011
#------------------------------------------------------------------------------

# Accepts a data.frame with columns CV, n, design  or CV, df
# if both n and df are present the df given take precedence
# returns a list with the components CV, df, CVupper, alpha
CVpooled <- function(CVdata, alpha=0.2, logscale=TRUE, robust=FALSE)
{
  if (nrow(CVdata)==0) {
    warning("No data in subset!")
    ret        <- list(CV=NA, df=0, CVupper=NA, alpha=alpha,
                       robust=robust)
    class(ret) <- "CVp"
    return(ret)
  }
  # some checking of input consistency
  if (!("CV" %in% names(CVdata))) stop("data.frame must have a column CV!")
  if (!("n" %in% names(CVdata)) & !("df" %in% names(CVdata))) {
    stop("Input data.frame must have a column n or df!")
  }
  if (!("design" %in% names(CVdata))) {
    CVdata$design <- "2x2"
    message("Classical 2x2 cross-over designs assumed for all entries.")
  }  
  if (!("df" %in% names(CVdata))) {
    CVdata$df <- 0
    for (i in seq_along(CVdata$n))
    { dno <- .design.no(CVdata$design[i])
      if (!is.na(dno)) {
        dprop <- .design.props(dno)
        n <- CVdata$n[i]
        if (robust) { 
          CVdata$df[i] <- eval(parse(text=dprop$df2,srcfile=NULL))
        } else {
          CVdata$df[i] <- eval(parse(text=dprop$df,srcfile=NULL))
        }
      } 
    }
  }
  # calculate se from CV
  CVdata$se  <- CVdata$CV
  if (logscale) CVdata$se  <- CV2se(CVdata$CV)
  # pooling of variance = se^2
  CVdata$se2 <- CVdata$se^2
  dftot      <- sum(CVdata$df, na.rm=TRUE)
  pooledse2  <- sum(CVdata$df*CVdata$se2, na.rm=TRUE)/dftot
  CVpooled   <- sqrt(pooledse2)
  if (logscale) CVpooled   <- se2CV(sqrt(pooledse2))
  # upper CL for CV
  chi        <- qchisq(alpha,dftot)
  CLCV       <- sqrt(pooledse2*dftot/chi)
  if (logscale) CLCV <- se2CV(CLCV)
  ret        <- list(CV=CVpooled, df=dftot, CVupper=CLCV, alpha=alpha, 
                     robust=robust)
  class(ret) <- "CVp"
  return(ret)
  
}

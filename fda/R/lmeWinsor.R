lmeWinsor <- function(fixed, data, random, 
        lower=NULL, upper=NULL, trim=0, quantileType=7,
        correlation, weights, subset, method, na.action,
        control, contrasts = NULL, keep.data=TRUE, ...)
{
##
## 1.  Identify inputs and outputs
##
  if(any(grep('lme4', search())>0))
    stop("lmeWinsor requires nlme, which can NOT be used",
         " with lme4 in the search path.") 
  library(nlme)
#
  cl <- match.call()
  if(missing(na.action))
    na.action <- get(options('na.action')$na.action) 
  mdly <- mdlx <- fixed
  mdly[[3]] <- NULL
  mdlx[[2]] <- NULL
#
  xNames <- unique(c(all.vars(mdlx), all.vars(random)))
  yName <- all.vars(mdly)
  if(length(yName) != 1)
    stop("'fixed' must include a single column of 'data'")
  if(as.character(mdly[[2]]) != yName)
    stop("lmWinsor can not accept a fixed with a transformed 'y';",
         "  left hand side = ", mdly[[2]], ";  y = ", yName)
##
## 2.  Check 'lower' and 'upper'
##
#  2.1.  Do lower and upper have names?    
  lowNames <- names(lower)
  if(length(lowNames) != length(lower))
    stop("lower must have names")
  hiNames <- names(upper)
  if(length(hiNames) != length(upper))
    stop("upper must have names") 
#  2.2.  Identify numeric columns of 'data' 
  numVars <- sapply(data, is.numeric)
  numV <- names(numVars[numVars])
#  2.3.  Are numeric variables in lower and upper?  
  numLower <- (numV %in% names(lower))
  numUpper <- (numV %in% names(upper))
  nnV <- length(numV)
#  2.4.  Some numeric variables are not in lower;  add   
  {
    if(!all(numLower)){
      Lower <- rep(NA, nnV)
      names(Lower) <- numV
      gotLo <- (numV %in% lowNames) 
      if(any(gotLo)) {
        loGot <- lower[numV[gotLo]]
        Lower[gotLo] <- loGot
      }
      for(v in numV[!numLower])
        Lower[v] <- quantile(data[[v]], trim, na.rm=TRUE, names=FALSE, 
                             type=quantileType)
    }
    else
      Lower <- lower
  }
  {
    if(!all(numUpper)){
      Upper <- rep(NA, nnV)
      names(Upper) <- numV
      gotHi <- (numV %in% hiNames)
      if(any(gotHi)){
        hiGot <- upper[numV[gotHi]]
        Upper[gotHi] <- hiGot
      }
      for(v in numV[!numUpper])
        Upper[v] <- quantile(data[[v]], 1-trim, na.rm=TRUE, names=FALSE, 
                             type=quantileType)
    }
    else
      Upper <- upper
  }
##
## 3.  clipData = data with xNames clipped to (Lower, Upper)
##
  clipData <- data
  for(x. in intersect(xNames, numV)){
    x.L <- Lower[x.]
    xl <- pmax(data[[x.]], x.L)
#    xl <- pmax(data[[x.]], x.L*(1+3*.Machine$double.eps))
    x.U <- Upper[x.]
    clipData[[x.]] <- pmin(xl, x.U) 
#    clipData[[x.]] <- pmin(xl, x.U *(1-3*.Machine$double.neg.esp))
  }
##
## 4.  fit <- lme(...)
##
  N <- nrow(data)
  if(missing(subset))subset <- 1:N 
#
  cl0 <- as.list(cl)
  cl0[[1]] <- NULL
  cl0$lower <- NULL
  cl0$upper <- NULL
  cl0$trim <- NULL
  cl0$quantileType <- NULL
  cl0$data <- as.name('clipData')
#
  fit <- do.call('lme', cl0)
##
## 5.  Convert to class 'lmWinsor'
##
  fit$call <- cl
  fit$lower <- Lower
  fit$upper <- Upper
  class(fit) <- c("lmeWinsor", class(fit))
##
## 6.  all fit$fitted %in% (Lower, Upper)[yName]?
##
  pred <- fit$fitted
  yL. <- Lower[yName]
#  
  yU. <- Upper[yName]
  pred[] <- pmax(yL., pmin(pred, yU.))
  fit$fitted <- pred
  resW <- (data[[yName]]-pred)
  deW <- (resW-fit$residuals) 
  fit$residuals <- resW
#
# fit$modelStruct
#
# fit$sigma 
#
##
## 7.  Done 
##
  fit
}

dtruncdist <- function(x, ..., dist='norm', truncmin=-Inf, 
                       truncmax=Inf){  
##
## 1.  dist
##
  ddist <- paste0('d', dist)
  dlist <- as.list(args(ddist))
  pdist <- paste0('p', dist) 
  plist <- as.list(args(pdist))
##  
## 2.  dots 
##
  dots <- list(...)
  pots <- dots
  pots$log <- NULL
##
## 3.  pdist(truncmin < X <= truncmax)
##
  if('log.p' %in% names(plist)){    
# Compute on the log scale by default 
# because that preserves numeric precision 
# for a broader range of values.
# If we keep unlogged values, we can get 
# underflow, the log of which would be -Inf.  
#    
#   log(x-n) = log(x) + log(1-n/x)
#     = log(x) + log(1-exp(log(n)-log(x)))    
    pots$log.p <- TRUE 
    pots$q <- truncmax
    lpmax <- do.call(pdist, pots)
    pots$q <- truncmin
    lpmin <- do.call(pdist, pots)
    lx_n <- (lpmax + log1p(-exp(lpmin-lpmax)))
#    x_n <-exp(lx_n)
  } else {
    pots$q <- truncmax 
    p.max <- do.call(pdist, pots)
    pots$q <- truncmin
    p.min <- do.call(pdist, pots)
    x_n <- (p.max-p.min)
    lx_n <- log(x_n)
  }
##
## 4.  ddist(x, ...)
##
  dots$x <- x 
  dDist <- do.call(ddist, dots) 
  if(('log' %in% names(dots)) && dots$log){
    out <- (dDist - lx_n)
#    out <- (dDist - log(x_n))
    out[x<=truncmin] <- (-Inf)
    out[truncmax<x] <- Inf 
  } else {
#    out <- (dDist / x_n) 
    out <- (dDist / exp(lx_n))
    out[x<=truncmin] <- 0 
    out[truncmax<x] <- 0 
  }
##
## 5.  done 
##
  out 
} 
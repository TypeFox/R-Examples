ptruncdist <- function(q, ..., dist='norm', truncmin=-Inf, 
                       truncmax=Inf){  
##
## 1.  dist
##
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
#   lx_n <- log(p.max - p.min)    
    lx_n <- (lpmax + log1p(-exp(lpmin-lpmax)))
#    x_n <-exp(lx_n)
  } else {
    pots$q <- truncmax 
    p.max <- do.call(pdist, pots)
    pots$q <- truncmin
    p.min <- do.call(pdist, pots)
    lpmin <- log(p.min)
    x_n <- (p.max-p.min)
    lx_n <- log(x_n)
  }
##
## 4.  pdist(q, ...)
##
  dots$q <- q 
  pDist <- do.call(pdist, dots) 
  if(('log' %in% names(dots)) && dots$log){
# log(pDi-x_n) = log(pDi) + log(1-x_n/pDi)
#   = log(pDi) + log1p(-exp(lx_n - pDist))  
    out <- rep_len(NA, length(pDist)) 
    out[q<=truncmin] <- (-Inf)
    out[truncmax<q] <- 0 
    pNo0 <- ((truncmin<q) & (q<=truncmax))
    ln.1mNn.Nz <- log1p(-exp(lpmin-pDist[pNo0]))
    out[pNo0] <- (pDist[pNo0] + ln.1mNn.Nz - lx_n)   
#    pD <- (pDist + log1p(-exp(lpmin-pDist)) - lx_n)
#    out <- (pD - lx_n)
#    out <- (dDist - log(x_n))
#    out[q<=truncmin] <- (-Inf)
#    out[truncmax<q] <- 0
  } else {
#    out <- (dDist / x_n) 
    x_n <- exp(lx_n)
    p.min <- exp(lpmin)
    out <- ((pDist-p.min) / x_n)
    out[q<=truncmin] <- 0 
    out[truncmax<q] <- 1
  }
##
## 5.  done 
##
  out 
} 
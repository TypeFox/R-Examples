qtruncdist <- function(p, ..., dist='norm', truncmin=-Inf, 
                       truncmax=Inf){  
##
## 1.  dist
##
  pdist <- paste0('p', dist)
  plist <- as.list(args(pdist))
  qdist <- paste0('q', dist) 
  qlist <- as.list(args(qdist))
##  
## 2.  dots 
##
  dots <- list(...)
  pots <- dots
  pots$log <- NULL
  if(('log.p' %in% names(dots)) && dots$log.p){ 
    oops <- (p>0)
    if(any(oops)){
      warning('log(p)>0 for ', sum(oops), 
            ' observations, the first of which is ', 
            which(oops)[1], '; will set to NaN')
    }
#   unlog p to simplify the algorithm  
#   Otherwise, the algorithm becomes so complicated
#     it's not worth the effort -- at least at the moment.
#   If we find an example where that complexity is needed, 
#     we can reconsider.
    p <- exp(p) # unlog 
    dots$log.p <- FALSE 
  } else {
    oops <- ((p<0) | (1<p))
    if(any(oops)){
      warning('probabilities outside [0, 1] for ', sum(oops), 
            ' observations, the first of which is ', 
            which(oops)[1], '; will set to NaN')
    }
  }
##
## 3.  pdist(truncmin < X <= truncmax)
##
  if('log.p' %in% names(plist)){
    pots$log.p <- FALSE
  }
  pots$q <- truncmax
  p.max <- do.call(pdist, pots)
  pots$q <- truncmin
  p.min <- do.call(pdist, pots)
  px_n <- (p.max-p.min)      
##
## 4.  qdist(p.min+px_n*p, ...)
##
  dots$p <- (p.min+px_n*p) 
  out <- do.call(qdist, dots) 
  out[oops] <- 0/0
##
## 5.  done 
##
  out 
} 
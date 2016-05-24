"isCOP.RTI" <-
function(cop=NULL, para=NULL, wrtV=FALSE, delta=0.005, ...) {
  T <- seq(0+delta, 1-delta, by=delta)
  is.RTI <- TRUE
  if(wrtV) {
    for(u in T) {
      derC  <- sapply(T, function(v) { return(derCOP2(u,v, cop=cop, para=para, ...))})
      CdivT <- sapply(T, function(v) { return(u - cop(u,v, para=para, ...)/(1 - v))})
      if(any(derC < CdivT)) { is.RTI <- FALSE; break }
    }
  } else {
    for(v in T) {
      derC  <- sapply(T, function(u) { return(derCOP(u,v, cop=cop, para=para, ...))})
      CdivT <- sapply(T, function(u) { return(v - cop(u,v, para=para, ...)/(1 - u))})
      if(any(derC < CdivT)) { is.RTI <- FALSE; break }
    }
  }
  return(is.RTI)
}

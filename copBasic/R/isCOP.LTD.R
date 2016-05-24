"isCOP.LTD" <-
function(cop=NULL, para=NULL, wrtV=FALSE, delta=0.005, ...) {
  T <- seq(0+delta, 1-delta, by=delta)
  is.LTD <- TRUE
  if(wrtV) {
    for(u in T) {
      derC  <- sapply(T, function(v) { return(derCOP2(u,v, cop=cop, para=para, ...))})
      CdivT <- sapply(T, function(v) { return(cop(u,v, para=para, ...)/v)})
      if(any(derC > CdivT)) { is.LTD <- FALSE; break }
    }
  } else {
    for(v in T) {
      derC  <- sapply(T, function(u) { return(derCOP(u,v, cop=cop, para=para, ...))})
      CdivT <- sapply(T, function(u) { return(cop(u,v, para=para, ...)/u)})
      if(any(derC > CdivT)) { is.LTD <- FALSE; break }
    }
  }
  return(is.LTD)
}

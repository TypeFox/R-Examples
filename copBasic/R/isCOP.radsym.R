"isCOP.radsym" <-
function(cop=NULL, para=NULL, delta=0.005, tol=1e-4, ...) {
  T <- seq(0+delta, 1-delta, by=delta)
  tmp <- sapply(T, function(u) { sapply(T, function(v) {
             return(u + v - 1 + cop(1-u,1-v, cop=cop, para=para, ...) -
                                cop(  u,  v, cop=cop, para=para, ...))
                               } ) })
  ifelse(any(abs(tmp) > tol), return(FALSE), return(TRUE))
}

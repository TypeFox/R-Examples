"COP" <-
function(u, v, cop=NULL, para=NULL, reflect=c("cop", "surv", "acute", "grave"), ...) {
   reflect <- match.arg(reflect)
   if(length(u) > 1 & length(v) > 1 & length(u) != length(v)) {
      warning("length u = ", length(u), " and length v = ", length(v))
      warning("longer object length is not a multiple of shorter object length, no recycling")
      return(NA)
   }
   if(length(u) == 1) {
      u <- rep(u, length(v))
   } else if(length(v) == 1) {
      v <- rep(v, length(u))
   }
   cop <- switch(reflect,
                 cop   =             cop(u,     v, para=para, ...),
                 surv  = u + v - 1 + cop(1-u, 1-v, para=para, ...),
                 acute =     v     - cop(1-u,   v, para=para, ...),
                 grave = u         - cop(  u, 1-v, para=para, ...) )
  return(cop)
}

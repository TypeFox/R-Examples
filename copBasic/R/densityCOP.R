"densityCOP" <-
function(u,v, cop=NULL, para=NULL,
              deluv=.Machine$double.eps^0.25, truncate.at.zero=TRUE, ...) {
   if (length(u) > 1 & length(v) > 1 & length(u) != length(v)) {
      warning("length u = ", length(u), " and length v = ", length(v))
      warning("longer object length is not a multiple of shorter object length, ",
              "no recycling")
      return(NA)
   }
   if(length(u) == 1) {
      u <- rep(u, length(v))
   } else if(length(v) == 1) {
      v <- rep(v, length(u))
   }
   den <- sapply(1:length(u), function(i) {
                     u1 <- u[i];     v1 <- v[i]
                     u2 <- u1+deluv; v2 <- v1+deluv
                     if(u2 > 1) { u1 <- u[i]-deluv; u2 <- u[i] }
                     if(v2 > 1) { v1 <- v[i]-deluv; v2 <- v[i] }
                     (cop(u2,v2, para=para, ...) -
                      cop(u2,v1, para=para, ...) -
                      cop(u1,v2, para=para, ...) +
                      cop(u1,v1, para=para, ...))/deluv^2 })
   if(truncate.at.zero) den[den < 0] <- 0
   return(den)
}


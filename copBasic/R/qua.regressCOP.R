"qua.regressCOP" <-
function(f=0.5, u=seq(0.01,0.99, by=0.01), cop=NULL, para=NULL, ...) {
  V <- sapply(1:length(u), function(i) {
          v <- derCOPinv(cop=cop, u[i], f, para=para, ...)
          if(is.na(v)) {
             warning("could not uniroot in derCOPinv, skipping i=",i," with U=", u[i])
             return(NA)
          }
          return(v) } )
  z <- data.frame(U=u,V=V)
  return(z)
}

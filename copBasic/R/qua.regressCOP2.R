"qua.regressCOP2" <-
function(f=0.5, v=seq(0.01, 0.99, by=0.01), cop=NULL, para=NULL, ...) {
  U <- sapply(1:length(v), function(i) {
          u <- derCOPinv2(cop=cop, v[i], f, para=para, ...)
          if(is.na(u)) {
             warning("could not uniroot in derCOPinv2, skipping i=", i," with V=", v[i])
             return(NA)
          }
          return(u) } )
  z <- data.frame(U=U,V=v)
  return(z)
}

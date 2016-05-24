"tailconCOP" <- function(t, cop=NULL, para=NULL, ...) {
   if(any(t < 0)) {
      warning("at least one t < 0, returning NULL")
      return(NULL)
   }
   if(any(t > 1)) {
      warning("at least one t > 1, returning NULL")
      return(NULL)
   }
   if (is.null(cop)) {
      warning("must have copula argument specified, returning NULL")
      return(NULL)
   }
   I <- ifelse(t < 0.5, 1, 0)
   delC <- COP(t, t, cop=cop, para=para, ...)
   ( delC/t * I) + ((1 - 2*t + delC)/(1-t) * ! I)
}

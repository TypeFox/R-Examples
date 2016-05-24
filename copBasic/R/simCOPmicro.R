"simCOPv" <-
function(u, cop=NULL, para=NULL, reflect=c("cop", "surv", "acute", "grave"), ...) {
   simCOPmicro(u, cop=cop, para=para, reflect=reflect, ...)
}


"simCOPmicro" <-
function(u, cop=NULL, para=NULL, reflect=c("cop", "surv", "acute", "grave"), ...) {
  reflect <- match.arg(reflect)
  n <- length(u); t <- runif(n); v <- vector(mode="numeric", length=n)
  v <- switch(reflect,
    cop   = sapply(1:n, function(i) {     derCOPinv(cop=cop,   u[i],   t[i], para=para, ...) }),
    surv  = sapply(1:n, function(i) { 1 - derCOPinv(cop=cop, 1-u[i], 1-t[i], para=para, ...) }),
    acute = sapply(1:n, function(i) {     derCOPinv(cop=cop, 1-u[i],   t[i], para=para, ...) }),
    grave = sapply(1:n, function(i) { 1 - derCOPinv(cop=cop,   u[i], 1-t[i], para=para, ...) })
  )
  if(any(is.na(v))) warning("could not uniroot at least for one element in derCOPinv")
  return(v)
}

# surv  is a reflection on the horizontal AND vertical axes
# acute is a reflection on the horizontal axis
# grave is a reflection on the verical axis

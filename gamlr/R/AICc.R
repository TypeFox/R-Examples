
## corrected AIC calculation
AICc <- function(object, k=2){
  ll <- logLik(object)
  d <- attributes(ll)$df
  n <- attributes(ll)$nobs
  ic <- -2*ll+ k*d*n/(n-d-1)
  ic[d+1 > n] <- Inf
  attributes(ic)[c('df','nobs','class')] <- NULL
  ic
}

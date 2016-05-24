extractAIC.manyglm<- function(fit, scale, k=2, ...)
{
  n.vars = NCOL(fit$y)
  ll=logLik(fit)
  edf = attr(ll,"df") * n.vars
  AIC = -2*sum(ll) + k * edf
  return( c(edf,AIC) )
  NextMethod("extractAIC")
}

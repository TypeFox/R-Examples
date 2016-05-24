AIC.evmOpt <- function(object, penalized=FALSE, ..., k=2){
  if (penalized) { ll <- object$ploglik }
  else {ll <- object$loglik }
  -2*ll + k*length(coef(object))
}

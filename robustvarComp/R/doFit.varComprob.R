doFit.varComprob  <- function(object) {
  if (isTRUE(object$doFit))
    return(object)
  ans <- eval(object$doFit)
  ans$na.action <- object$na.action
  ans$offset <- ans$offset
  ans$contrasts <- object$contrasts
  ans$xzlevels <- object$xzlevels
  ans$call <- object$call
  ans$terms <- object$terms
  ans$model <- object$model
  ans$X <- object$X
  ans$Y <- object$Y
  ans$K <- object$K
  ans$random.labels <- object$random.labels
  ans$weights <- object$weights
  return(ans)	
}

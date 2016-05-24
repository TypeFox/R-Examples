coef.mcemGLMM <- function(object,...) {
  coef0 <- tail(object$mcemEST, n = 1)[1:ncol(object$x)]
  names(coef0) <- colnames(object$mcemEST)[1:ncol(object$x)]
  return(coef0)
}
ranef.mcemGLMM <- function(object, ...) {
  return(colMeans(object$randeff))
}
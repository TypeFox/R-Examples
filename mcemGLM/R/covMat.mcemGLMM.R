covMat.mcemGLMM <- function(object, ...) {
  kP <- ncol(object$x)
  fMat <- solve(object$iMatrix)
  fMat <- fMat[1:kP, 1:kP]
  colnames(fMat) <- colnames(object$x)
  rownames(fMat) <- colnames(object$x)
  return(fMat)
}
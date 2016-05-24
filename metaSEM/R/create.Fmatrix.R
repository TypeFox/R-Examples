create.Fmatrix <- function(x, name, as.mxMatrix=TRUE, ...) {
  x <- as.logical(x)
  Fmatrix <- Diag(as.numeric(x))[x, , drop=FALSE]
  if (as.mxMatrix) {
    if (missing(name)) as.mxMatrix(Fmatrix) else as.mxMatrix(Fmatrix, name=name)
  } else {
    Fmatrix
  }
}

matrix2bdiag <- function(x, ...) {
  tmp <- split(as.matrix(x), row(x))
  # Use bdiagMat() to handle string matrices
  out <- bdiagMat(lapply(tmp, vec2symMat, ...))
  as.matrix(out)
}

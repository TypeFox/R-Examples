is.diag <- function (c, EPS=1e-12) {
  #  tests a matrix for being diagonal
  #  C   ... matrix to be tested
  #  EPS ... testing criterion: max off-diagonal element over min diagonal
  #          element must be less than EPS
  if (!is.matrix(c)) return(FALSE)
  cd <- dim(c)
  if (cd[1] != cd[2]) return(FALSE)
  mindg <- min(abs(diag(c)))
  maxodg <- max(abs(c - diag(diag(c))))
  if (maxodg/mindg < EPS) return(TRUE) else return(FALSE)
}

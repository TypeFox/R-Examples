`outer.2s` <-
function (x, y, f, dn=TRUE, ...) {
  if (!(is.vector(x) && is.vector(y) && is.numeric(x) && is.numeric(y))) {
    stop("arguments not numeric vectors")
  }
  result <- matrix(NA, nrow=length(x), ncol=length(y))
  for (i in 1:length(x)) {
    for (j in 1:length(y)) {
      result[i,j] <- f(x[i], y[j], ...)
      }
    }
  if (dn) {
    lab.x                   <- deparse(substitute(x))
    lab.y                   <- deparse(substitute(y))
    dimnames(result)        <- list(x, y)
    names(dimnames(result)) <- c(lab.x, lab.y)
  }
  return(result)
}


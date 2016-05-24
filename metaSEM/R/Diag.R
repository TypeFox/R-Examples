Diag <- function(x, ...) {
#Diag <- function(x, nrow, ncol) {
  if (class(x)=="character" & missing(...)) {
    p <- length(x)
    out <- matrix(0, nrow=p, ncol=p)
    diag(out) <- x
  } else {
    out <- diag(x, ...)
    #out <- diag(x, nrow=nrow, ncol=ncol)
  }
  out
}

`Diag<-` <- function(x, value) {
  diag(x) <- value
  x
}

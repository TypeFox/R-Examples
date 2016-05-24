
rvmatrix <- function (data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL) {
  data <- as.vector(data)
  if (missing(nrow)) 
      nrow <- ceiling(length(data)/ncol)
  else if (missing(ncol)) 
      ncol <- ceiling(length(data)/nrow)
  if (byrow) {
    X <- t(rvarray(data, c(ncol, nrow), dimnames=dimnames))
  } else {
    X <- rvarray(data, c(nrow, ncol), dimnames=dimnames)
  }
  return(X)
}

rvarray <- function (data = NA, dim = length(data), dimnames = NULL) {
  as.rv(array(data = data, dim = dim, dimnames = dimnames))
}

is.matrix.rv <- function (x) {
  dx <- dim(x)
  return((!is.null(dx)) && length(dx)==2)
}

as.matrix.rv <- function (x, ...) {
  if (is.matrix(x)) {
    return(x)
  }
  dn <- if (!is.null(names(x))) list(names(x), NULL) else NULL
  rvarray(x, dim=c(length(x), 1), dimnames=dn)
}


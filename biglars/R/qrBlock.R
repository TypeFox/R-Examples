
# 'R' and 'Qty' behave as in/out variables.
#

qrBlock <- function(X, y, R, Qty)
{
  nrowX <- nrow(X)
  ncolX <- ncol(X)
  dimR <- ncol(R)

  Ry <- .Fortran("rblock",
      as.integer(dimR),
      R = R,
      as.integer(dimR),
      as.integer(nrowX),
      as.integer(ncolX),
      X,
      as.integer(nrowX),
      double(nrowX),
      double(dimR),
      Qty = as.double(Qty),
      as.double(y),
      PACKAGE = "biglars")[c('R', 'Qty')]
}


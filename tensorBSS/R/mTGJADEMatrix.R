mTGJADEMatrix <- function(x, i, j, lags){
  RES <- .Call( "mTGJADEMatrixC", x, i, j, lags, PACKAGE = "tensorBSS")
  RES
}

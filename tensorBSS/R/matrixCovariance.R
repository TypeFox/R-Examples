matrixCovariance <- function(x){
  RES <- .Call( "matrixCovarianceC", x, PACKAGE = "tensorBSS")
  RES
}

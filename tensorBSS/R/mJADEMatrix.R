mJADEMatrix <- function(x, i, j, cov){
  RES <- .Call( "mJADEMatrixC", x, i, j, cov, PACKAGE = "tensorBSS")
  RES
}

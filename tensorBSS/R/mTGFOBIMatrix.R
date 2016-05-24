mTGFOBIMatrix <- function(x, lag){
  RES <- .Call( "mTGFOBIMatrixC", x, lag, PACKAGE = "tensorBSS")
  RES
}

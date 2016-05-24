mFOBIMatrix <- function(x){
  RES <- .Call( "mFOBIMatrixC", x, PACKAGE = "tensorBSS")
  RES
}

mAutoCovMatrix <- function(x, lag){
  RES <- .Call( "mAutoCovMatrixC", x, lag, PACKAGE = "tensorBSS")
  RES
}

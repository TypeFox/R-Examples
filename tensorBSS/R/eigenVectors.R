eigenVectors <- function(x){
  RES <- .Call( "eigenVectorsC", x, PACKAGE = "tensorBSS")
  RES
}

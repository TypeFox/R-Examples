symmetricPower <- function(x, r){
  RES <- .Call( "symmetricPowerC", x, r, PACKAGE = "tensorBSS")
  RES
}

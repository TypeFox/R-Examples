"singcheck.bigqr" <-
function(bigQR){
  bigQR <- .Call("singcheckQR", bigQR)
}


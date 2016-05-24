varlz <-
function(vz,z){
# varlz() - variance of log of z
# vlz <- log(1 + vz/(z*z))
  tol <- 1.e-8
  tol <- 1 + tol
  cvp <- 1 + vz/(z*z)
  if(cvp < tol) {
#   vlz <- log(tol)
    vlz <- cvp - 1
  }
  else{
    vlz <- log(cvp)
  }
  return(vlz)
}

covlyz <-
function(cyz, y, z){
# covlyz() - covariance of log y and log z
# clyz <- log(1 + cyz/(y*z))
  tol <- 1.e-8
  ccp <- 1 + cyz/(y*z)
  if(ccp < tol) {
#   clyz <- log(tol)
    clyz <- ccp -1
  }
  else {
    clyz <- log(ccp)
  }
  return(clyz)
}

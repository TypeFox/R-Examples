varz <-
function(vlz,z){
# varz() - variance of z given variance of logz
# vz <- z * z *(exp(vlz) - 1)
  p <- exp(vlz)
  if(!is.na(p)){
    if(p < 1) {
      vz <- z * z * p
    }
    else {
      vz <- z * z * (p - 1)
    }
  }
  else {
    vz <- NA
  }
  return(vz)
}

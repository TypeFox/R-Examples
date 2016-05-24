
ldist=function( z, q0 = q0, lambda0) {
  term1 = (1.+z)^2
  term2 =  1.+2.*(q0+lambda0)*z
  term3 = z*(2.+z)*lambda0
  denom = (term1*term2 - term3)
  out = z*0.
  good =which(denom>0.0)
  out[good] = 1./sqrt(denom[good])
  return(out)
}
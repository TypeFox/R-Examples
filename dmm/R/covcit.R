covcit <-
function(v,l,iv, vsiga, il,jl){
#  covcit()  - sampling covariance of component iv var with total var
  covcitot <- 0
  ijb <- (il-1)*l+jl
  for(ic in 1 : v) {
    covcitot <- covcitot + vsiga[(ijb-1)*v+ic,(ijb-1)*v+iv]  # sums only row iv (vx1)
  }
  return(covcitot)
}

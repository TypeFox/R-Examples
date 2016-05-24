vt <-
function(v,l,vsiga, il, jl){
# vt(() - sampling variance of total of v components for traits il,jl
  vtot <- 0
  ijb <- (il-1)*l+jl
  for(iv in 1 : v) {
    for(jv in 1:v) {
      vtot <- vtot + vsiga[(ijb-1)*v+iv,(ijb-1)*v+jv]  # sums whole block(vxv)
    }
  }
  return(vtot)
}

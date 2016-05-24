covt <-
function(v,l,vsiga,il,jl,kl){
# covt()  - sampling cov of phenotypic var &/or cov's
  covtot <- 0
  ijb <- (il-1)*l+jl
  kb <- (kl-1)*l+kl
  for(iv in 1:v){
    for(jv in 1:v){
      covtot <- covtot + vsiga[(ijb-1)*v+iv, (kb-1)*v+jv] # sums whole block for traits ij and k
    }
  }
  return(covtot)
}

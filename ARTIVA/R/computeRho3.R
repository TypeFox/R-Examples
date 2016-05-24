computeRho3 <-
function(k, kmin, kmax, c, lambda){
  # INPUT:l= k the number of hidden states
  #          kmax the maximal number of hidden states
  #          c constante
  #          lambda simulation parameter for the numbr of hidden states
  # OUTPUT: rho
  # depends on: (/)

  rho=array(1,3)
  if(k == kmax) { rho[1]=0 } else { rho[1] = c*min(1,lambda/(k+1)) }
  if(k == kmin) { rho[2] = rho[1] } else { rho[2] = rho[1]+c*min(1,k/lambda) }

  return(rho)
}

betaIntEst <-
function(locs, var, L0, L){
  
  lambda <- {(locs - L0) * (L - locs)} / var  - 1
  
  ### Estimating alpha  
  alpha <- lambda * {(locs - L0) / (L - L0)}
  
  ### Estimating beta
  beta <- lambda * {(L - locs) / (L - L0)}
  
  return (list(alpha.init = alpha, beta.init = beta))  
}

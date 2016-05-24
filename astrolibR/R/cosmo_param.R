cosmo_param = function(omega_m, omega_lambda, omega_k, q0) {
  
  nk = missing(omega_k)
  nlambda = missing(omega_lambda)
  nomega = missing(omega_m)
  nq0 = missing(q0)
  if((nk + nlambda + nomega + nq0)>2 )
    stop('at least 2 cosmological parameters must be specified')

  
  if((!nomega && !nlambda) ){ 
    if(nk) omega_k = 1 - omega_m - omega_lambda 
    if(nq0) q0 = omega_m/2. - omega_lambda
  }
  else if((!nomega && !nk) ){ 
    if(nlambda) omega_lambda = 1. -omega_m - omega_k
    if(nq0) q0 = -1 + omega_k + 3*omega_m/2
  }
  else if((!nlambda && !nk) ){ 
    if(nomega) omega_m = 1.-omega_lambda - omega_k
    if(nq0) q0 = (1 - omega_k - 3.*omega_lambda)/2.
  }
  else if((!nomega && !nq0) ){
    if(nk) omega_k = 1 + q0 - 3*omega_m/2. 
    if(nlambda) omega_lambda  = 1. - omega_m - omega_k
  }
  else if((!nlambda && !nq0) ){
    if(nk) omega_k = 1 - 2*q0 - 3*omega_lambda
    if(nomega) omega_m = 1.-omega_lambda - omega_k
  }
  else if((!nk && !nq0) ){
    if(nomega) omega_m = (1 + q0 - omega_k)*2/3.
    if(nlambda) omega_lambda = 1. - omega_m - omega_k
  }

  ##?? arnab: how can length(omega_k)==0 now??
  
  if(length(omega_k)==0 ) omega_k = 0 #default is flat space
  if(length(omega_lambda)==0 ) omega_lambda = 0.7
  if(length(omega_m)==0 ) omega_m = 1 - omega_lambda
  if(length(q0)==0 ) q0 = (1 - omega_k - 3*omega_lambda)/2.
  return(list(omega_m=omega_m, omega_lambda=omega_lambda,
              omega_k=omega_k, q0=q0))
}

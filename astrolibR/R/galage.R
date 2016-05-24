#source('cosmo_param.R')

galage = function (z, zform, h0=70.0, omega_m, lambda0, k,q0, silent=FALSE) {

  tmp = cosmo_param(omega_m, lambda0, k, q0)
  omega_m = tmp$omega_m
  omega_lambda = tmp$omega_lambda
  k = tmp$omega_k
  q0 = tmp$q0
  
  nz = length(z)
  age = z*0.            #return same dimensions and data type as z
                                        #
  
  for(i in 1:nz) {
    if((z[i]>=zform) )
      age_z = 0
    else 
    
      age_z = integrate('dtdz', z[i], zform, q0 = q0, lambda0 =
            lambda0)$value
    cat('age_z = ',age_z)
    age[i] = age_z
  }
  return (age * 3.085678e+19 / 3.15567e+7/ h0)
}


lumdist=function(z, h0=70, k=0, lambda0=0.7, omega_m=0.3, q0=0.55) {

  n = length(z)
  tmp = cosmo_param(omega_m,lambda0, k, q0)
  omega_m = tmp$omega_m
  lambda0 = tmp$omega_lambda
  k = tmp$omega_k
  q0 = tmp$q0

  c = 2.99792458e5                  #  speed of light in km/s

  if(missing(h0)) h0=70

  if(lambda0== 0 ){
    denom = sqrt(1+2*q0*z) + 1 + q0*z
    dlum = (c*z/h0)*(1 + z*(1-q0)/denom)
    return(dlum)
  }

    dlum = z*0.0
    for( i in 1:n) {
      if(z[i]<=0.0 ) {
        dlum[i] = 0.0
      }
      else {
        lz = integrate(ldist,0,z[i],q0 = q0, lambda0 = lambda0)$value
        dlum[i] = lz
      }
    }
    if(k>0 )
      dlum = sinh(sqrt(k)*dlum)/sqrt(k)
    else if(k<0 )
      dlum = sin(sqrt(-k)*dlum)/sqrt(-k) > 0

    return(c*(1+z)*dlum/h0)

}

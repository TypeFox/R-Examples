jprecess = function(ra, dec,  mu_radec,
  parallax=numeric(n), rad_vel=numeric(n), epoch=1950) {

  n = length(ra)
  if(length( rad_vel )!=n ) {
    stop(paste('rad_vel keyword vector must contain',n,'values'))
  }
  
  if(!missing(mu_radec)){
    if((length( mu_radec)!=2*n ))
    stop(paste('mu_radec keyword (proper motion) be dimensioned (2',n,')'))
  }

  radeg = 180/pi
  sec_to_radian = 1./radeg/3600.0e0
  m =  cbind( c(+0.9999256782, +0.0111820610, +0.0048579479,  
          -0.000551,     +0.238514,     -0.435623     ), 
        c( -0.0111820611, +0.9999374784, -0.0000271474,   
         -0.238565,     -0.002667,      +0.012254      ), 
        c( -0.0048579477, -0.0000271765, +0.9999881997 , 
         +0.435739,      -0.008541,      +0.002117      ), 
        c( +0.00000242395018, +0.00000002710663, +0.00000001177656, 
         +0.99994704,    +0.01118251,    +0.00485767    ), 
        c( -0.00000002710663, +0.00000242397878, -0.00000000006582, 
         -0.01118251,     +0.99995883,    -0.00002714    ), 
        c( -0.00000001177656, -0.00000000006587, 0.00000242410173, 
         -0.00485767,   -0.00002718,     1.00000956) )
  a = 1e-6*c(-1.62557, -0.31919, -0.13843)        #in radians
  a_dot = 1e-3*c(1.244, -1.579, -0.660)           #in arc seconds per century
  if(epoch!=1950.0 )
    a = a + sec_to_radian * a_dot * (epoch - 1950.0)/100.0
  ra_rad = ra/radeg       ;      dec_rad = dec/radeg
  cosra =  cos( ra_rad )  ;       sinra = sin( ra_rad )
  cosdec = cos( dec_rad ) ;      sindec = sin( dec_rad )
  ra_2000 = ra*0.
  dec_2000 = dec*0.
  for(i  in 1:n) {
    r0 = c( cosra[i]*cosdec[i], sinra[i]*cosdec[i], sindec[i] )
    if(missing( mu_radec) ){
      mu_a = 0.0e0
      mu_d = 0.0e0
    } 
    else {
      mu_a = mu_radec[ (2*n-1)]
      mu_d = mu_radec[ 2*n]
    }
    r0_dot = c( -mu_a*sinra[i]*cosdec[i] - mu_d*cosra[i]*sindec[i],  #velocity vector
              mu_a*cosra[i]*cosdec[i] - mu_d*sinra[i]*sindec[i] , 
              mu_d*cosdec[i] )  + 21.095 * rad_vel[i] * parallax[i] * r0
                                        # remove the effects of the e-terms of aberration to form r1 and r1_dot.
    r1 = r0 - a + (sum(r0 * a))*r0
    r1_dot = r0_dot - a_dot + ( sum( r0 * a_dot))*r0
    r_1 = c(r1, r1_dot)         
    
    r = m %*% r_1
    if(missing(mu_radec) ){
      rr = r[1:3]
      v =  r[4:6]
      t = ((epoch - 1950.0e0) - 50.00021)/100.0e0
      rr1 = rr + sec_to_radian*v*t
      x = rr1[1]  
      y = rr1[2]  
      z = rr1[3]  
    } else {
      x = r[1]  
      y = r[2]  
      z = r[3]  
      x_dot = r[4]   
      y_dot = r[5]  
      z_dot = r[6]
    }
    r2 = x^2 + y^2 + z^2
    rmag = sqrt( r2 )
    dec_2000[i] = asin( z / rmag)
    ra_2000[i] = atan2( y, x)
    if(!missing(mu_radec) ){
      mu_radec[ (2*n-1) ] = ( x*y_dot - y*x_dot) / ( x^2 + y^2)
      mu_radec[ 2*n ] = ( z_dot* (x^2 + y^2) - z*(x*x_dot + y*y_dot) ) /  
        ( r2*sqrt( x^2 + y^2) )
    }
    
    if(parallax[i]>0 ){
      rad_vel[i] = ( x*x_dot + y*y_dot + z*z_dot )/ (21.095*parallax[i]*rmag)
      parallax[i] = parallax[i] / rmag
    }
  } 
  neg = (ra_2000<0)
  ra_2000[neg] = ra_2000[neg] + 2*pi
  ra_2000 = ra_2000*radeg
  dec_2000 = dec_2000*radeg

  return(list(ra_2000 = ra_2000, dec_2000 = dec_2000))
}

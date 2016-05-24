bprecess = function(
  ra,
  dec,
  mu_radec, 
  parallax = numeric(length(ra)),
  rad_vel = numeric(length(ra)),
  epoch = 2000) {

  n = length( ra ) 

  if(length(rad_vel)!=n )
        stop(paste('rad_vel keyword vector must contain ',
                    n,' values'))
 
 if(!missing(mu_radec) && (length(mu_radec)!=2*n ))
   stop('mu_radec keyword (proper motion) be dimensioned (2,' + 
                  n, ')')
  
  radeg = 180/pi
  sec_to_radian = 1/radeg/3600
 m =  cbind( c(+0.9999256795, -0.0111814828, -0.0048590040,  
         -0.000551,  -0.238560,     +0.435730)     , 
        c(+0.0111814828, +0.9999374849, -0.0000271557,   
         +0.238509,     -0.002667,      -0.008541)     , 
        c(+0.0048590039, -0.0000271771, +0.9999881946 , 
         -0.435614,      +0.012254,      +0.002117)      , 
        c(-0.00000242389840, +0.00000002710544, +0.00000001177742, 
         +0.99990432,    -0.01118145,    -0.00485852)    , 
        c(-0.00000002710544, -0.00000242392702, +0.00000000006585, 
         +0.01118145,     +0.99991613,    -0.00002716)    , 
        c(-0.00000001177742, +0.00000000006585,-0.00000242404995, 
         +0.00485852,   -0.00002717,    +0.99996684)) 
 a_dot = 1e-3*c(1.244, -1.579, -0.660 )           #in arc seconds per century
  ra_rad = ra/radeg
  dec_rad = dec/radeg
  cosra =  cos( ra_rad )
  sinra = sin( ra_rad )
  cosdec = cos( dec_rad )
  sindec = sin( dec_rad )
  dec_1950 = dec*0.
  ra_1950 = ra*0.
  for(i  in 1:n) {
    a = 1e-6*c( -1.62557, -0.31919, -0.13843)        #in radians
    r0 = c( cosra[i]*cosdec[i], sinra[i]*cosdec[i], sindec[i] )
    if(!missing(mu_radec) ){
      mu_a = mu_radec[ (2*n-1) ]
      mu_d = mu_radec[ 2*n ]
      r0_dot = c( -mu_a*sinra[i]*cosdec[i] -
                mu_d*cosra[i]*sindec[i] ,  #velocity vector
                mu_a*cosra[i]*cosdec[i] -
                mu_d*sinra[i]*sindec[i] , 
                mu_d*cosdec[i] ) +
                  21.095 * rad_vel[i] * parallax[i] * r0
    }
    else {
      r0_dot = c(0.0, 0.0, 0.0)
    }

    r_0 = c(r0, r0_dot)
    
    r_1 =  r_0 %*% t(m)
                                        # include the effects of the e-terms of aberration to form r and r_dot.
    r1 = r_1[1:3]  
    r1_dot = r_1[4:6] 
    if(!missing(mu_radec) ){
      r1 = r1 + sec_to_radian * r1_dot * (epoch - 1950.0)/100.
      a = a + sec_to_radian * a_dot * (epoch - 1950.0)/100.
 }
    x1 = r_1[1]
    y1 = r_1[2]
    z1 = r_1[3]
    rmag = sqrt( x1^2 + y1^2 + z1^2 )
    s1 = r1/rmag
    s1_dot = r1_dot/rmag
    s = s1

    for(j  in 1:3) {
      r = s1 + a - (sum(s * a))*s
      s = r/rmag
    } 
    x = r[1]
    y = r[2]
    z = r[3]  
    r2 = x^2 + y^2 + z^2
    rmag = sqrt( r2 )
 
    if(!missing(mu_radec) ){
         r_dot = s1_dot + a_dot - ( sum( s * a_dot))*s
         x_dot = r_dot[1]  ; y_dot= r_dot[2]  ;  z_dot = r_dot[3]
         mu_radec[(2*n-1)] = ( x*y_dot - y*x_dot) / ( x^2 + y^2)
         mu_radec[2*n] = ( z_dot* (x^2 + y^2) - z*(x*x_dot + y*y_dot) ) /  
                     ( r2*sqrt( x^2 + y^2) )
    }
    dec_1950[i] = asin( z / rmag)
    ra_1950[i] = atan2( y, x)
    if(parallax[i]>0 ){
      rad_vel[i] = ( x*x_dot + y*y_dot + z*z_dot )/ (21.095*parallax[i]*rmag)
      parallax[i] = parallax[i] / rmag
    }
 }
 neg = (ra_1950<0)
 ra_1950[neg] = ra_1950[neg] + 2.*pi
 ra_1950 = ra_1950*radeg 
 dec_1950 = dec_1950*radeg

 return(list(ra_1950 = ra_1950, dec_1950=dec_1950))
}

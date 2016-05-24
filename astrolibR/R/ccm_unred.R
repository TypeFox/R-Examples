ccm_unred = function( wave, flux, ebv, R_V=3.1) {

  x = 10000./ wave                # Convert to inverse microns 
  npts = length( x )
  a = numeric(npts)  
  b = numeric(npts)
  good =which( (x>0.3) & (x <1.1)) #Infrared
  a[good] =  0.574 * x[good]^(1.61)
  b[good] = -0.527 * x[good]^(1.61)

  good =which( (x>=1.1) & (x<3.3) )
                                        #Use new constants from O'Donnell (1994)
  y = x[good] - 1.82
  c1 = c( 1. , 0.104,   -0.609,    0.701,  1.137,        #New coefficients
    -1.718,   -0.827,    1.647, -0.505 )        #from O'Donnell
  c2 = c( 0.,  1.952,    2.908,   -3.989, -7.985,        #(1994)
    11.102,    5.491,  -10.805,  3.347 )
  a[good] = polyidl( y, c1)
  b[good] = polyidl( y, c2)
  
  good =which( (x>=3.3) & (x<8) );Ngood=length( good )  #Mid-UV
  if(Ngood>0 ){
    y = x[good]
    F_a = numeric(Ngood)
    F_b = numeric(Ngood)
    good1 =which( (y>5.9))

    y1 = y[good1] - 5.9
    F_a[ good1] = -0.04473 * y1^2 - 0.009779 * y1^3
    F_b[ good1] =   0.2130 * y1^2  +  0.1207 * y1^3
    
    
    a[good] =  1.752 - 0.316*y - (0.104 / ( (y-4.67)^2 + 0.341 )) + F_a
    b[good] = -3.090 + 1.825*y + (1.206 / ( (y-4.62)^2 + 0.263 )) + F_b
  }
  
  good =which( (x>=8) & (x<=11)) #Far-UV

  y = x[good] - 8.
  c1 = c( -1.073, -0.628,  0.137, -0.070 )
  c2 = c( 13.670,  4.257, -0.420,  0.374 )
  a[good] = polyidl(y, c1)
  b[good] = polyidl(y, c2)
  
  A_V = R_V * ebv
  A_lambda = A_V * (a + b/R_V)
  funred = flux * 10.^(0.4*A_lambda)       #Derive unreddened flux
  return(funred)     
}                               

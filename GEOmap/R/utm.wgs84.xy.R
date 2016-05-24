`utm.wgs84.xy` <-
  function(  phideg,  lamdeg, PROJ.DATA )
  {
############
    
#########   latD=14; longD=-91;  long0D=264
#########    http://www.uwgb.edu/dutchs/UsefulData/UTMFormulas.htm


      long0D = PROJ.DATA$LON0

      latD = phideg
      longD =     lamdeg
      
      lat =   latD*pi/180
#### long =  RPMG::fmod(longD, 360)*pi/180
#### long0 = RPMG::fmod(long0D, 360)*pi/180
      
      
      long =  longD*pi/180
      long0 = long0D*pi/180
      
      
      a=6378137  	
    b=6356752.3142  
    f = (a-b)/a	
###f=1/298.257223563



    k0 = 0.9996 
    ##   = scale along long0 
    e = sqrt(1-b^2/a^2) 

    ## .08 approximately. This is the eccentricity of the earth's elliptical cross-section.

    ep2 = (e*a/b)^2

###  e^2/(1-e^2) 
    ##   = .007 approximately
###  The quantity e' only occurs in even powers so it need only be calculated as e'2.

    n = (a-b)/(a+b)
    rho = a*(1-e^2)/((1-(e*sin(lat))^2)^(3/2))   
    ##  . This is the radius of curvature of the earth in the meridian plane.
    
    nu = a/sqrt(1-(e*sin(lat))^2)

                                        # . This is the radius of curvature of the earth perpendicular to the meridian plane. It is also the distance 
                                        #          	from the point in question to the polar axis, measured perpendicular to the earth's surface.

    p = (long-long0)



    Ap = a*(1 - n + (5/4)*(n^2 - n^3) + (81/64)*(n^4 - n^5) )

    Bp = (3*a*n/2)*(1 - n + (7/8)*(n^2 - n^3) + (55/64)*(n^4 - n^5))

    Cp = (15*a*n^2/16)*(1 - n + (3/4)*(n^2 - n^3) )

    Dp = (35*a*n^3/48)*(1 - n + (11/16)*(n^2 - n^3) )

    Ep = (315*a*n^4/51)*(1 - n )

    S = Ap*lat - Bp*sin(2*lat) + Cp*sin(4*lat) - Dp*sin(6*lat) + Ep*sin(8*lat) 


########## , where lat is in radians and 
    

    M = a *( (1 - e^2/4 - 3*e^4/64 - 5*e^6/256 )*lat 
      - (3*e^2/8 + 3*e^4/32 + 45*e^6/1024)*sin(2*lat) 
      + (15*e^4/256 + 45*e^6/1024 )*sin(4*lat)
      - (35*e^6/3072 )* sin(6*lat) )



    K1 = S*k0

    K2 = k0 *nu *sin(lat)*cos(lat)/2

    K3 = (k0*nu *sin(lat)*((cos(lat))^3) /24) * (5 - (tan(lat))^2 + 9*ep2*(cos(lat))^2 + 4*ep2^2*(cos(lat))^4)
    

    K4 = k0*nu *cos(lat)

    K5 = (k0*nu *((cos(lat))^3)/6)*(1 - (tan(lat))^2 + ep2*(cos(lat))^2)


    y =  K1 + K2*p*p + K3*p^4

    x = 500000+  K4*p + K5*p^3          ############  , where


    return(list(x=x, y=y))


  }

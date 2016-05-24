`UTM.ll` <-
  function( x, y, PROJ.DATA)
  {
### Converting UTM to Latitude and Longitude

### y = northing, x = easting
###     (relative to central meridian; subtract 500,000 from conventional UTM coordinate).
### Calculate the Meridional Arc

### This is easy:

    long0 = PROJ.DATA$LON0*pi/180

    xprime = x-500000
    if(is.na(PROJ.DATA$DATUM$a))
      {
        a=6378137
        b=6356752.3142
      }
    else
        {
          a=PROJ.DATA$DATUM$a
          b=PROJ.DATA$DATUM$b
        }
    
    
    f = (a-b)/a	
###f=1/298.257223563

    k0 = 0.9996 
    M = y/k0
### Calculate Footprint Latitude
    e = sqrt(1-b^2/a^2) 

    esqred = e^2

    mu = M/( a*(1 - esqred/4 - 3*e^4/64 - 5*e^6/256)  )
    e1 = (1 - sqrt(1 - esqred) )/(1 + sqrt(1 - esqred))

###footprint latitude

###   , where:

    J1 = (3*e1/2 - 27*e1^3/32 )
    J2 = (21*e1^2/16 - 55*e1^4/32 )
    J3 = (151*e1^3/96 )
    J4 = (1097*e1^4/512 )

    
    fp = mu + J1*sin(2*mu) + J2*sin(4*mu) + J3*sin(6*mu) + J4*sin(8*mu)
###Calculate Latitude and Longitude

    ep2 = (e*a/b)^2
                              ##      = esqred/(1-esqred) 
    C1 = ep2*(cos(fp)^2)
    T1 = (tan(fp))^2
    R1 = a*(1-esqred)/(1-esqred* (sin(fp))^2  )^(3/2) 

###This is the same as rho in the forward conversion formulas above,
####  but calculated for fp instead of lat.

    N1 = a/(1-esqred*(sin(fp))^2)^(1/2)

###This is the same as nu in the forward conversion formulas above,
    ##   but calculated for fp instead of lat.

    D = xprime/(N1*k0)


###, where:

    Q1 = N1 *tan(fp)/R1
    Q2 = (D^2/2)
    Q3 = (5 + 3*T1 + 10*C1 - 4*C1^2 -9*ep2)*D^4/24
    Q4 = (61 + 90*T1 + 298*C1 +45*T1^2  - 3*C1^2 -252*ep2)*D^6/720
    lat = fp - Q1*(Q2 - Q3 + Q4)


###, where:

    Q5 = D
    Q6 = (1 + 2*T1 + C1)*D^3/6
    Q7 = (5 - 2*C1 + 28*T1 - 3*C1^2 + 8*ep2 + 24*T1^2)*D^5/120

    long = long0 + (Q5 - Q6 + Q7)/cos(fp)

    lon = long*180/pi
    lat =  lat*180/pi

    return(list(lat=lat, lon=lon))


    


  }

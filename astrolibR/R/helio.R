#source('cirrange.R')
helio=function(jd, list1, radian=FALSE) {

  pd = cbind( c( 0.38709893, 0.20563069, 7.00487,  48.33167,  77.45645, 252.25084 ), 
        c( 0.72333199, 0.00677323, 3.39471,  76.68069, 131.53298, 181.97973 ),  
        c( 1.00000011, 0.01671022, 0.00005, -11.26064, 102.94719, 100.46435), 
        c( 1.52366231, 0.09341233, 1.85061,  49.57854, 336.04084, 355.45332), 
        c( 5.20336301, 0.04839266, 1.30530, 100.55615,  14.75385,  34.40438),  
        c( 9.53707032, 0.05415060, 2.48446, 113.71504,  92.43194,  49.94432), 
        c(19.19126393, 0.04716771, 0.76986,  74.22988, 170.96424, 313.23218),  
        c(30.06896348, 0.00858587, 1.76917, 131.72169,  44.97135, 304.88003), 
        c(39.48168677, 0.24880766,17.14175, 110.30347,
        224.06676, 238.92881))
  
  dpd = cbind(  c(0.00000066, 0.00002527, -23.51, -446.30, 573.57, 538101628.29 ), 
         c( 0.00000092, -0.00004938, -2.86, -996.89, -108.80, 210664136.06), 
         c(-0.00000005, -0.00003804, -46.94, -18228.25, 1198.28, 129597740.63),  
         c(-0.00007221, 0.00011902, -25.47, -1020.19, 1560.78, 68905103.78 ), 
         c(0.00060737, -0.00012880, -4.15, 1217.17, 839.93, 10925078.35 ), 
         c(-0.00301530, -0.00036762, 6.11, -1591.05, -1948.89, 4401052.95),  
         c(0.00152025, -0.00019150, -2.09, -1681.40, 1312.56, 1542547.79 ), 
         c(-0.00125196, 0.0000251, -3.64, -151.25, -844.43, 786449.21 ), 
         c(-0.00076912, 0.00006465, 11.07, -37.33, -132.25, 522747.90) ) 
  jd0 = 2451545.0    #julian date for epoch 2000.0
  radeg = 180/pi
  t = (jd - jd0)/36525.0          #time in centuries since 2000.0
  
  ip = list1
  dpd[3:6,ip] = dpd[3:6,ip]/3600.0       #convert arc seconds to degrees
  ntime = length(t)
  nplanet = length(list1)
  hrad = matrix(0,nplanet,ntime)
  hlong = hrad
  hlat = hrad
  
  for(i  in 1:ntime) {         #sml made longword
    pd1 = pd[,ip,drop=F] + dpd[,ip,drop=F]*t[i]
    a = pd1[1,]                            #semi-major axis
    eccen = pd1[2,]                        #eccentricity
    n = 0.9856076686/a/sqrt(a)/radeg      #mean motion, in radians/day
    l =  pd1[6,]/radeg                     #mean longitude
    ppi = pd1[5,]/radeg                  #longitude of the perihelion
    omega = pd1[4,]/radeg               #longitude of the ascending node
    inc = pd1[3,]/radeg          #inclination in radians
   

    m = l - ppi
    m=cirrange(m,radians=T)

    e1 = m + (m + eccen*sin(m) - m)/(1 - eccen*cos(m) )
    e = e1 + (m + eccen*sin(e1) - e1)/(1 - eccen*cos(e1) )
    maxdif = max(abs(e-e1))
    niter = 0
    while ((maxdif>1e-5) && (niter<10)) {        
      e1 = e
      e = e1 + (m + eccen*sin(e1) - e1)/(1 - eccen*cos(e1) )
      maxdif = max(abs(e-e1))
      niter = niter+1
    }
    
    
    nu = 2*atan( sqrt( (1+eccen)/(1-eccen) )* tan(e/2))   #true anomaly
    hrad[,i] =  as.numeric(a*(1 - eccen*cos(e) )) 
    hlong[,i] = as.numeric (nu + ppi)                 
    hlat[,i] = as.numeric( asin(sin(hlong[,i] - omega)*sin(inc) ) )
  } 
  hlong = cirrange(hlong,radians=T) 
  if(!radian){ 
    hlong = hlong*radeg
    hlat = hlat*radeg
  }

  return(list(hrad=hrad,hlong=hlong, hlat=hlat))  
}

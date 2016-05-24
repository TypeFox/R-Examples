#source('nutate.r')
#source('cirrange.r')
#source('ten.r')
#source('polyidl.r')

moonpos=function( jd, radian=F){

  npts = length(jd)
  dtor = pi/180.0
                                        #  form time in julian centuries from 1900.0
  t = (jd - 2451545.0)/36525.0e0
  d_lng = c(0,2,2,0,0,0,2,2,2,2,0,1,0,2,0,0,4,0,4,2,2,1,1,2,2,4,2,0,2,2,1,2,0,0, 
    2,2,2,4,0,3,2,4,0,2,2,2,4,0,4,1,2,0,1,3,4,2,0,1,2,2)
  m_lng = c(0,0,0,0,1,0,0,-1,0,-1,1,0,1,0,0,0,0,0,0,1,1,0,1,-1,0,0,0,1,0,-1,0, 
    -2,1,2,-2,0,0,-1,0,0,1,-1,2,2,1,-1,0,0,-1,0,1,0,1,0,0,-1,2,1,0,0)
  mp_lng = c(1,-1,0,2,0,0,-2,-1,1,0,-1,0,1,0,1,1,-1,3,-2,-1,0,-1,0,1,2,0,-3,-2,
    -1,-2,1,0,2,0,-1,1,0,-1,2,-1,1,-2,-1,-1,-2,0,1,4,0,-2,0,2,1,-2,-3,2,1,-1, 
    3,-1)
  f_lng = c(0,0,0,0,0,2,0,0,0,0,0,0,0,-2,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0, 
    0,0,0,-2,2,0,2,0,0,0,0,0,0,-2,0,0,0,0,-2,-2,0,0,0,0,0,0,0,-2)
  sin_lng = c(6288774,1274027,658314,213618,-185116,-114332,58793,57066,53322, 
    45758,-40923,-34720,-30383,15327,-12528,10980,10675,10034,8548,-7888,-6766, 
    -5163,4987,4036,3994,3861,3665,-2689,-2602,2390,-2348,2236,-2120,-2069,2048, 
    -1773,-1595,1215,-1110,-892,-810,759,-713,-700,691,596,549,537,520,-487, 
    -399,-381,351,-340,330,327,-323,299,294,0.0)
  cos_lng = c(-20905355,-3699111,-2955968,-569925,48888,-3149,246158,-152138, 
    -170733,-204586,-129620,108743,104755,10321,0,79661,-34782,-23210,-21636, 
    24208,30824,-8379,-16675,-12831,-10445,-11650,14403,-7003,0,10056,6322, 
    -9884,5751,0,-4950,4130,0,-3958,0,3258,2616,-1897,-2117,2354,0,0,-1423, 
    -1117,-1571,-1739,0,-4421,0,0,0,0,1165,0,0,8752.0)
  d_lat = c(0,0,0,2,2,2,2,0,2,0,2,2,2,2,2,2,2,0,4,0,0,0,1,0,0,0,1,0,4,4,0,4,2,2,
    2,2,0,2,2,2,2,4,2,2,0,2,1,1,0,2,1,2,0,4,4,1,4,1,4,2)
  m_lat = c(0,0,0,0,0,0,0,0,0,0,-1,0,0,1,-1,-1,-1,1,0,1,0,1,0,1,1,1,0,0,0,0,0,0,
    0,0,-1,0,0,0,0,1,1,0,-1,-2,0,1,1,1,1,1,0,-1,1,0,-1,0,0,0,-1,-2)
  mp_lat = c(0,1,1,0,-1,-1,0,2,1,2,0,-2,1,0,-1,0,-1,-1,-1,0,0,-1,0,1,1,0,0,3,0,  
    -1,1, -2,0,2,1,-2,3,2,-3,-1,0,0,1,0,1,1,0,0,-2,-1,1,-2,2,-2,-1,1,1,-1,0,0)
  f_lat =c( 1,1,-1,-1,1,-1,1,1,-1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,3,1,1,1,-1,-1,-1, 
    1,-1,1,-3,1,-3,-1,-1,1,-1,1,-1,1,1,1,1,-1,3,-1,-1,1,-1,-1,1,-1,1,-1,-1, 
    -1,-1,-1,-1,1)
  sin_lat = c(5128122,280602,277693,173237,55413,46271,32573,17198,9266,8822, 
    8216,4324,4200,-3359,2463,2211,2065,-1870,1828,-1794,-1749,-1565,-1491, 
    -1475,-1410,-1344,-1335,1107,1021,833,777,671,607,596,491,-451,439,422, 
    421,-366,-351,331,315,302,-283,-229,223,223,-220,-220,-185,181,-177,176, 
    166,-164,132,-119,115,107.0)
  coeff0 = c(218.3164477, 481267.88123421, -0.0015786e0, 1.0/538841.0, 
    -1.0/6.5194e7 )
  lprimed = polyidl(t, coeff0)
  lprimed = cirrange(lprimed)
  lprime = lprimed*dtor
  coeff1 = c(297.8501921, 445267.1114034, -0.0018819, 1.0/545868.0, 
    -1.0/1.13065e8 )
  d = polyidl(t, coeff1)
  d = cirrange(d)
  d = d*dtor
  coeff2 = c(357.5291092, 35999.0502909, -0.0001536, 1.0/2.449e7 )
  m = polyidl(t,coeff2) 
  m=cirrange(m) 
  m = m*dtor
  coeff3 = c(134.9633964, 477198.8675055, 0.0087414, 1.0/6.9699e4, 
    -1.0/1.4712e7 )
  mprime = polyidl(t, coeff3) 
  mprime = cirrange(mprime)
  mprime = mprime*dtor
  coeff4 = c(93.2720950, 483202.0175233, -0.0036539, -1.0/3.526e7, 
    1.0/8.6331e8 )
  f = polyidl(t, coeff4 ) 
  f = cirrange(f)
  f = f*dtor
  e = 1 - 0.002516*t - 7.4e-6*t^2
  e2 = e^2
  ecorr1 = which(abs(m_lng)==1)
  ecorr2 = which(abs(m_lat)==1)
  ecorr3 = which(abs(m_lng)==2)
  ecorr4 = which(abs(m_lat)==2)
  a1 = (119.75 + 131.849*t) * dtor
  a2 = (53.09 + 479264.290*t) * dtor
  a3 = (313.45 + 481266.484*t) * dtor
  suml_add = 3958*sin(a1) + 1962*sin(lprime - f) + 318*sin(a2)
  sumb_add =  -2235*sin(lprime) + 382*sin(a3) + 175*sin(a1-f) +  
    175*sin(a1 + f) + 127*sin(lprime - mprime) - 
      115*sin(lprime + mprime)
  geolong = numeric(npts)
  geolat = geolong
  dis = geolong
  for(i in 1:npts) {
    sinlng = sin_lng ; coslng = cos_lng ; sinlat = sin_lat
    sinlng[ecorr1] = e[i]*sinlng[ecorr1]
    coslng[ecorr1] = e[i]*coslng[ecorr1]
    sinlat[ecorr2] = e[i]*sinlat[ecorr2]
    sinlng[ecorr3] = e2[i]*sinlng[ecorr3]
    coslng[ecorr3] = e2[i]*coslng[ecorr3]
    sinlat[ecorr4] = e2[i]*sinlat[ecorr4]
    arg = d_lng*d[i] + m_lng*m[i] +mp_lng*mprime[i] + f_lng*f[i]
    geolong[i] = lprimed[i] + ( sum(sinlng*sin(arg)) + suml_add[i] )/1.0e6
    dis[i] = 385000.56 + sum(coslng*cos(arg))/1.0e3
    arg = d_lat*d[i] + m_lat*m[i] +mp_lat*mprime[i] + f_lat*f[i]
    geolat[i] = (sum(sinlat*sin(arg)) + sumb_add[i])/1.0e6
    
  }
  tmp = nutate(jd)
  nlong = tmp$nut_long
  elong = tmp$nut_obliq                     #find the nutation in longitude
  geolong= geolong + nlong/3.6e3
  geolong = cirrange(geolong)
  lambda = geolong*dtor
  beta = geolat*dtor
  c = c(21.448,-4680.93,-1.55,1999.25,-51.38,-249.67,-39.05,7.12,27.87,5.79,2.45)
  epsilon = ten(23,26) + polyidl(t/1.e2,c)/3600.
  eps = (epsilon + elong/3600. )*dtor          #true obliquity in radians
  #cat('epsilon=',epsilon,'\n')
  #cat('elong=',elong,'\n')
  #cat('dtor=',dtor,'\n')
  ra = atan2( sin(lambda)*cos(eps) - tan(beta)* sin(eps), cos(lambda) )
  #cat('lambda=',lambda,'\n')
  #cat('eps=',eps,'\n')
  #cat('beta=',lambda,'\n')

  ra = cirrange(ra,radians=T)
  dec = asin( sin(beta)*cos(eps) + cos(beta)*sin(eps)*sin(lambda) )

  if(!radian){
    ra = ra/dtor
    dec = dec/dtor
  } else {
    geolong = lambda
    geolat = beta
  }
return(list(ra=ra, dec=dec, dis=dis, geolong=geolong, geolat=geolat))
}

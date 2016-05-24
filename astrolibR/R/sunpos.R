sunpos=function( jd, radian=F) {

  dtor = pi/180.0       #(degrees to radian, as.eouble precision)
  t = (jd - 2415020.0)/36525.0e0
  L = (279.696678e0+((36000.768925e0*t) %% 360.0e0))*3600.0e0
  me = 358.475844e0 + ((35999.049750e0*t) %% 360.0e0)
  ellcor  = (6910.1e0 - 17.2e0*t)*sin(me*dtor) + 72.3e0*sin(2.0e0*me*dtor)
  L = L + ellcor
  mv = 212.603219e0 + ((58517.803875e0*t) %% 360.0e0) 
  vencorr =
    4.8e0 * cos((299.1017e0 + mv - me)*dtor) + 
      5.5*cos((148.3133e0 + 2.0e0 * mv - 2.0e0 * me )*dtor) + 
      2.5*cos((315.9433e0 + 2.0e0 * mv - 3.0e0 * me )*dtor) + 
      1.6*cos((345.2533e0 + 3.0e0 * mv - 4.0e0 * me )*dtor) + 
      1.0*cos((318.15e0   + 3.0e0 * mv - 5.0e0 * me )*dtor)
  
  L = L + vencorr
  mm = 319.529425e0  +  (( 19139.858500e0 * t) %%  360.0)
  marscorr =
    2*cos((343.8883e0 -2*mm  +  2*me)*dtor) + 
    1.8*cos((200.4017e0 -  2*mm  + me)*dtor)
  
  L = L + marscorr
  mj = 225.328328e0  +  (( 3034.6920239e0 * t)  %%  360.0e0 )
  jupcorr = 7.2e0 * cos(( 179.5317e0 - mj + me )*dtor) + 
    2.6e0 * cos((263.2167e0  -  mj ) *dtor) + 
      2.7e0 * cos(( 87.1450e0  -  2.0e0 * mj  +  2.0e0 * me ) *dtor) + 
        1.6e0 * cos((109.4933e0  -  2.0e0 * mj  +  me ) *dtor)
  L = L + jupcorr
  d = 350.7376814e0  + (( 445267.11422e0 * t)  %%  360.0e0 )
  mooncorr  = 6.5e0 * sin(d*dtor)
  L = L + mooncorr
  longterm  = + 6.4e0 * sin(( 231.19e0  +  20.20e0 * t )*dtor)
  L  =    L + longterm
  L  =  ( L + 2592000.0e0)  %%  1296000
  longmed = L/3600.0e0
  L  =  L - 20.5e0
  omega = 259.183275e0 - (( 1934.142008e0 * t ) %% 360.0e0 )
  L  =  L - 17.2e0 * sin(omega*dtor)
  oblt  = 23.452294e0 - 0.0130125e0*t + (9.2e0*cos(omega*dtor))/3600
  L = L/3600
  ra  = atan2( sin(L*dtor) * cos(oblt*dtor) , cos(L*dtor) )

  neg =which(ra<0.0e0) 
  ra[neg] = ra[neg] + 2.0*pi
  dec = asin(sin(L*dtor) * sin(oblt*dtor))
  
  if(radian){
    oblt = oblt*dtor 
    longmed = longmed*dtor
  } else {
    ra = ra/dtor
    dec = dec/dtor
  }

  return(list(ra=ra,dec=dec,longmed=longmed,oblt=oblt))
}

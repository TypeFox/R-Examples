`AlongGreat` <-
function( LON1,  LAT1,  km1, ang,  EARTHRAD = 6371 )
{

                                        #     /*
                                        #        given a point (lat lon) degrees
                                        #        a distance in km (km1) and
                                        #       an azimuthal direction Az
                                        #        calculate the point (phi, lam) which lies angle radians away from (lat, lon)
                                        #     */
                                        #

if(missing(EARTHRAD)){ EARTHRAD = 6371 }
  
  DEG2RAD = pi/180;
  RAD2DEG = 180/pi;
  phi1= DEG2RAD*LAT1;
  lam0= DEG2RAD* LON1;
  Az = DEG2RAD*ang;
  c = km1/EARTHRAD;
  
  phi =  asin( sin(phi1)*cos(c) + cos(phi1)*sin(c)*cos(Az));
  tem = atan2( sin(c)*sin(Az), (cos(phi1)*cos(c) - sin(phi1)*sin(c)*cos(Az)));
  lam = lam0 + tem;
  ddeg = RAD2DEG*c;

  lam = RAD2DEG*lam
  phi = RAD2DEG*phi
  return(list(lat=phi, lon=lam , distdeg=ddeg , distkm=c));
}


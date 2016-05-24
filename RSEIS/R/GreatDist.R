`GreatDist` <-
function( LON1,  LAT1,  LON2,  LAT2,  EARTHRAD = 6371)
{
 #  /* calculate the great arc distance 
 #     between two points
 #     given the lat lons in degrees */
 #
  if(missing(EARTHRAD)){ EARTHRAD = 6371 }
   DEG2RAD = pi/180;
   RAD2DEG = 180/pi;
  phi1= DEG2RAD*LAT1;
  lam0=DEG2RAD* LON1;
  phi=DEG2RAD*LAT2;
  lam=DEG2RAD*LON2;
  
  s1 = sin( (phi-phi1 )/2);
  s2 = sin( (lam-lam0 )/2);
  tem = sqrt(s1*s1 + cos(phi1)*cos(phi)*s2*s2);
  tem2 = 2*asin(tem);
  # tem2 is now in radians, convert to degrees and km
  
  tem3  =  RAD2DEG*tem2;
  tem4 = tem2*EARTHRAD;
  
  return(list(drad=tem2, ddeg=tem3, dkm=tem4));
}


`to.spherical` <-
function(x, y, z)
  {
###########  this matches C-code used in SDRSUBS.c
##########  /* Convert from Cartesian to spherical (az, dip in degrees) coordinates */
    RAD2DEG = 180/pi;
     rlength = sqrt( x*x + y*y + z*z);
     if(any(rlength==0.0)) return(list(az=NaN,dip=NaN, x=NaN, y=NaN, z=NaN));
     diprad = asin(z/rlength);
     azrad = atan2(y, x);
     az = azrad * RAD2DEG;
     dip = diprad * RAD2DEG;
     
     return(list(az=az,dip=dip, x=x, y=y, z=z));
   }


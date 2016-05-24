`TOSPHERE.DIP` <-
function(x, y, z)
  {
#####  /* Convert from Cartesian to spherical (az, dip in degrees) coordinates */
    RAD2DEG = 180/pi;
  length = sqrt( x*x + y*y + z*z);
  if(any(length==0.0)) return(list(az=NaN,dip=NaN, x=NaN, y=NaN, z=NaN));
  diprad = asin(z/length);
  azrad = atan2(y, x);
  az = azrad * RAD2DEG;
  dip = diprad * RAD2DEG;
  if(dip<0)
    {
      dip = 90-dip
      az = az+180
    }
  return(list(az=az,dip=dip, x=x, y=y, z=z))
}


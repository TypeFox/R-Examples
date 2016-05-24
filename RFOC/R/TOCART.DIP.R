`TOCART.DIP` <-
function(az, dip)
{
#####   Convert from spherical (az, dip in degrees) to Cartesian coordinates
#####     x pos north, y pos east, z pos downward
#####    az in degrees from north, dip in degrees from Z (down)
  DEG2RAD = pi/180;
  azrad = az * DEG2RAD;
  diprad = dip * DEG2RAD;
  
  z = sin(diprad);
  temp = cos(diprad);
  x = cos(azrad) * temp;
  y = sin(azrad) * temp;
  len=sqrt(x*x+y*y+z*z)
  z = z/len
  x = x/len
  y = y/len 
  return(list(x=x,y=y,z=z, az=az, dip=dip));
}


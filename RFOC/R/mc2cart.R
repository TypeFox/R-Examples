`mc2cart` <-
function(az, dip)
{
   DEG2RAD = pi/180;
    azrad = az * DEG2RAD;
    diprad = dip * DEG2RAD;
    z = sin (diprad);
    temp = cos(diprad);
    x = cos(azrad) * temp;
    y = sin(azrad) * temp;
    return(list(x=x, y=y, z=z))
}


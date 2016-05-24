`tosphereL` <-
function(A)
{
#######  /* Convert from Cartesian to spherical (az, dip in degrees) coordinates */
  a = TOSPHERE(A$x, A$y, A$z)
  return(a)
}


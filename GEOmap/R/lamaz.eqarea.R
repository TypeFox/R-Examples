`lamaz.eqarea` <-
function( phi1,  lam0,  phi,  lam)
{
  R= 1
  
  kprime = sqrt(2/(1+sin(phi1)*sin(phi)+cos(phi1)*cos(phi)*cos(lam-lam0) ) );
  
  x = R*kprime*cos(phi)*sin(lam-lam0);
  y = R*kprime*(cos(phi1)*sin(phi)-sin(phi1)*cos(phi)*cos(lam-lam0));
  
  return(list(x=x, y=y))
}


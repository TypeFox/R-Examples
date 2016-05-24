`ternfoc.point` <-
function(deltaB,  deltaP,  deltaT)
{

######### /* plot a focal point on a ternary diagram  */
  DEG2RAD = pi/180


  db=deltaB*DEG2RAD;
  dp=deltaP*DEG2RAD;
  dt=deltaT*DEG2RAD;
  
  censin=sin(DEG2RAD*35.26);
  cencos=cos(DEG2RAD*35.26);
  
  phi =atan2( sin(dt),sin(dp)) - 45*DEG2RAD;

  denom = (censin*sin(db)+cencos*cos(db)*cos(phi) );
  
  h = (cos(db)*sin(phi))/denom;
  
  v = (cencos*sin(db)-censin*cos(db)*cos(phi))/denom;

  return(list(h=h, v=v))
}


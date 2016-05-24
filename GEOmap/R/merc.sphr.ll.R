`merc.sphr.ll` <-
function(LON0, x, y)
  {
    R.MAPK= MAPconstants()$R.MAPK
    RAD2DEG=180/pi
    
    
    phi = 90 - RAD2DEG * 2 * atan(exp(-y / R.MAPK));
    lam = RAD2DEG * (x / R.MAPK) + LON0;
    return(list(lat=phi, lon=lam))
  }


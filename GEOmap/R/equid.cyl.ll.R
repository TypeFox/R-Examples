`equid.cyl.ll` <-
function(LON0, phi0, x, y)
  {

   R.MAPK= MAPconstants()$R.MAPK
   DEG2RAD=pi/180
   RAD2DEG=180/pi
   
    phi = RAD2DEG * y/R.MAPK
    lam =RAD2DEG * (x / (R.MAPK*cos(DEG2RAD *phi0))) + LON0;
    return(list(lat=phi, lon=lam))
  }


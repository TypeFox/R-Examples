`utm.sphr.ll` <-
function( x, y, PROJ.DATA)
  {
    k0 = 1.0;
    R<-	6378.2064
    RAD2DEG =180/pi
    DEG2RAD=pi/180

    R.MAPK=MAPconstants()$R.MAPK
    
    x = x+ PROJ.DATA$FE
    y = y+ PROJ.DATA$FN

    
    D = y/(R.MAPK*k0) + DEG2RAD *PROJ.DATA$LAT0;
    a1 = RAD2DEG *asin( sin(D)/cosh(x/(R.MAPK*k0)) );
    
    a2 = RAD2DEG * atan( sinh(x/(R.MAPK*k0))/cos(D)) + PROJ.DATA$LON0;
    
    phi = a1;
    lam = a2;
    
    return(list(lat=a1, lon=a2))
  }


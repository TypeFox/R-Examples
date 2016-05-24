`utm.sphr.xy` <-
function(  phi,  lam, PROJ.DATA )
  {
   ############  snyder  page 58
      DEG2RAD=pi/180
      R.MAPK = MAPconstants()$R.MAPK
      
    k0 = 1.0;
    theta = (lam - PROJ.DATA$LON0)
    
    B = cos(DEG2RAD *phi) * sin(DEG2RAD *theta);
    x = 0.5*R.MAPK*k0*log( (1+B)/(1-B));

    a = tan(DEG2RAD *phi)
    b = cos(DEG2RAD *theta)
    
    aa = atan( a/b)
    ww = which(b<0)

      si = sign(a[ww])
      si[si==0] = 1
    aa[ww] = aa[ww]+si*pi

   ## aa[phi== 0.0 &   ww ] = aa[phi== 0.0 & ww]+ pi
    
    y = R.MAPK * k0* ( aa  - DEG2RAD*PROJ.DATA$LAT0 );

    x = x - PROJ.DATA$FE
    y = y - PROJ.DATA$FN

    
    return(list(x=x, y=y))
  }


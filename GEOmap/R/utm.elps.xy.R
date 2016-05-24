`utm.elps.xy` <-
function(  phideg,  lamdeg, PROJ.DATA )
  {
   ############  snyder  page 58
    
    k0 = 0.9996;
    e2 = MAPconstants()$E2.MAPK
    ep2 = e2/(1-e2)
    
    DEG2RAD=pi/180
    RAD2DEG=180/pi

    A.MAPK=MAPconstants()$A.MAPK
    
    theta = DEG2RAD *(lamdeg - PROJ.DATA$LON0)
    phi = DEG2RAD *phideg
    
    N = A.MAPK/sqrt((1-e2*sin(phi)^2))
    Tee = tan(phi)^2
    C  = ep2*(cos(phi)^2)
    A  = theta*cos(phi)

    
 ###   M  = A.MAPK*( (1-e2/4 - 3*e2^2/64 -5*e2^3/256)*phi-(3*e2/8+3*e2^2/32
 ###     +45*e2^3/1024)*sin(2*phi)+(15*e2^2/256+45*e2^3/1024)*sin(4*phi)-(35*e2^3/3072)*sin(6*phi))

    ###   for clark can use simplified formula:
    M = 111132.0894*phideg-16216.94*sin(2*phi)+17.21*sin(4*phi)-0.02*sin(6*phi)

    M0 = 111132.0894*PROJ.DATA$LAT0-16216.94*sin(2*DEG2RAD*PROJ.DATA$LAT0)+17.21*sin(4*DEG2RAD*PROJ.DATA$LAT0)-0.02*sin(6*DEG2RAD*PROJ.DATA$LAT0)
 
    x = k0*N*(A+(1-Tee+C)*A^3/6+(5-18*Tee+Tee^2+72*C-58*ep2)*A^5/120)
    y = k0*(M-M0+N*tan(phi)*(A^2/2+(5-Tee+9*C+4*C^2)*A^4/24+(61-58*Tee+Tee^2+600*C-330*ep2)*A^6/720))

     x = x + PROJ.DATA$FE 
     y = y + PROJ.DATA$FN
    
   
    return(list(x=x, y=y))
  }


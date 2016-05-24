`utm.elps.ll` <-
function( x, y, PROJ.DATA)
  {

    k0 = 0.9996;
    e2 = MAPconstants()$E2.MAPK
    ep2 = e2/(1-e2)

     DEG2RAD=pi/180
    RAD2DEG=180/pi
    
    A.MAPK=MAPconstants()$A.MAPK
    

    phi0 = DEG2RAD*PROJ.DATA$LAT0
    lam0 = DEG2RAD*PROJ.DATA$LON0

    
    x = x-PROJ.DATA$FE
    y = y-PROJ.DATA$FN
 
    M0 = 111132.0894*PROJ.DATA$LAT0-16216.94*sin(2*DEG2RAD*PROJ.DATA$LAT0)+17.21*sin(4*DEG2RAD*PROJ.DATA$LAT0)-0.02*sin(6*DEG2RAD*PROJ.DATA$LAT0)
    M = M0 +y/k0
    
    e1 = (1-sqrt(1-e2))/(1+sqrt(1-e2))
    mu = M/(A.MAPK*(1-e2/4-3*e2^2/64-5*e2^3/256))
    
    phi1 = mu + (3*e1/2-27*e1^3/32)*sin(2*mu)+(21*e1^2/16-
      55*e1^4/32)*sin(4*mu)+(151*e1^3/96)*sin(6*mu)+(1097*e1^4/512)*sin(8*mu)
    
   C1 = ep2*cos(phi1)^2
    T1 = tan(phi1)^2
    N1 = A.MAPK/sqrt(1-e2*sin(phi1)^2)
    R1 = A.MAPK*(1-e2)/(1-e2*sin(phi1)^2)^(1.5)
    D = x/(N1*k0)

     phi = phi1 -(N1*tan(phi1)/R1)*(D^2/2-(5+3*T1+10*C1-4*C1^2-9*ep2)*D^4/24
       + (61+90*T1+298*C1+45*T1^2 - 252*ep2-3*C1^2)*D^6/720)

       lam = lam0 +(D - (1+2*T1 +C1)*D^3/6 + (5-2*C1+28*T1
          -3*C1^2 +8*ep2+24*T1^2)*D^5/120)/cos(phi1)

    a1 = RAD2DEG*phi
    a2 = RAD2DEG*lam
   
    return(list(lat=a1, lon=a2))
  }


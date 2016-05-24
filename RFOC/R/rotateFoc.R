 TP2XYZ<-function(trend,  dip)
  {
    t = trend * 0.017453293;
    p = dip * 0.017453293;

    x = c(cos(t) * cos(p),
      sin(t) * cos(p),
      sin(p))
    return(x)
  }

RotTP<-function(rotmat, strk1, dip1)
  {
    

    v1 = TP2XYZ(strk1, dip1)
    
    v2 = rotmat %*% v1;
    
### if the vector is pointing down, reflect it through zero
###     this is to get the upper hemisphere projection of poles */

    if(v2[3]< 0.0)
      {
        v2 = -1*v2
     }

    dip1 = abs((asin(v2[3])))

    if(v2[2]  != 0.0) {  strk1 = atan2(v2[2], v2[1]) }
    else { strk1 = pi/2 }

    dip1 =  57.295778667*dip1;
    strk1 =  57.295778667*strk1;

    return(list(strk=strk1, dip=dip1))
  }


JMAT<-function(phi)
  {
####  phi in degrees
    m1 = ROTY(90)
    m2 = ROTX(phi)
    JJ = m2 %*% m1
    JJ = JJ[1:3, 1:3]
    return( JJ )
  }

Rotfocphi<-function(phi, urot, udip, vrot,
                    vdip, az1, d1, az2, d2,
                    prot, pdip, trot, tdip)
  {
    

    p = phi;

    
    rotmat =  JMAT(p);

   ### print(rotmat)
    
### /* Dump_Mat(rotmat); */

    A1 = RotTP(rotmat, urot, udip);

    urot=A1$strk
    udip=A1$dip

    
    A2 = RotTP(rotmat, vrot, vdip);
    vrot=A2$strk
    vdip=A2$dip


    az1 =  urot - 180.0;
    d1 = 90.0 - udip;
    
    az2 =  vrot + 180.0;
    d2 = 90.0 - vdip;

    if(az1 < 0.0) az1 = az1 + 360.0;
    if(az2 < 0.0) az2 =  az2 + 360.0;

    P =   RotTP(rotmat, prot, pdip);

    prot = P$strk
    pdip = P$dip
    
    TEE =   RotTP(rotmat, trot, tdip);
    trot = TEE$strk
    tdip = TEE$dip
    
####    GG = c(urot, udip, vrot,
####      vdip, az1, d1, az2, d2,
####      prot, pdip, trot, tdip)

####

        GG = list(U=list(az=urot, dip=udip), V=list(az= vrot,
      dip=vdip), A1 =list(az=az1, dip=d1), A2=list(az=az2, dip=d2),
     P=list(az=prot, dip=pdip) , T=list(az=trot, dip=tdip))

    return(GG)
  }















rotateFoc<-function(MEX, phi)
  {
    vertrot<-function(phi)
      {
####  phi in degrees
        m1 = ROTY(-90)
        m2 = ROTX(-phi)
        return( m2 %*% m1 )
      }

    Rmat =  vertrot(phi) 

    P =     TOCART.DIP(MEX$P$az, MEX$P$dip)
    T=      TOCART.DIP(MEX$T$az, MEX$T$dip)
    U =     TOCART.DIP(MEX$U$az, MEX$U$dip)
    V =     TOCART.DIP(MEX$V$az, MEX$V$dip)
    F =     TOCART.DIP(MEX$F$az, MEX$F$dip)
    G =     TOCART.DIP(MEX$G$az, MEX$G$dip)

    
    A1 = TOCART.DIP(MEX$az1, MEX$dip1)
    A2 = TOCART.DIP(MEX$az2, MEX$dip2)

    vA1 = c(A1$x, A1$y, A1$z, 1)    
    rA1 = Rmat %*% vA1

    vA2 = c(A2$x, A2$y, A2$z, 1)    
    rA2 = Rmat %*% vA2

    vP = c(P$x, P$y, P$z, 1)    
    rP = Rmat %*% vP

    vT = c(T$x, T$y, T$z, 1)    
    rT = Rmat %*% vT 

    vU = c(U$x, U$y, U$z, 1)    
    rU = Rmat %*% vU

    vV = c(V$x, V$y, V$z, 1)    
    rV = Rmat %*% vV

    vF = c(F$x, F$y, F$z, 1)    
    rF = Rmat %*% vF

    vG = c(G$x, G$y, G$z, 1)    
    rG = Rmat %*% vG


    P  = TOSPHERE.DIP(rP[1,1], rP[2,1],rP[3,1] )
    T = TOSPHERE.DIP(rT[1,1], rT[2,1],rT[3,1] )
    V = TOSPHERE.DIP(rV[1,1], rV[2,1],rV[3,1] )
    U = TOSPHERE.DIP(rU[1,1], rU[2,1],rU[3,1] )
    G = TOSPHERE.DIP(rG[1,1], rG[2,1],rG[3,1] )
    F = TOSPHERE.DIP(rF[1,1], rF[2,1],rF[3,1] )
    
    
    ang2 = GetRakeSense(U$az, U$dip, V$az, V$dip, P$az, P$dip,  T$az, T$dip)

    M1= GetRake(F$az-90, F$dip,  G$az -90,  G$dip, ang2)

    mc = CONVERTSDR(M1$az1, M1$dip1, M1$rake1 )



    
    MEC = MRake(mc$M)
    MEC$UP = FALSE
    MEC$icol =  foc.icolor(MEC$rake1)
    MEC$ileg =  focleg(MEC$icol)
    MEC$fcol =   foc.color(MEC$icol)
    MEC$CNVRG = NA
    MEC$LIM = c(-1,-1, +1, +1)

    return(MEC)

  }

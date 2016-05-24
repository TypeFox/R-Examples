`Bfocvec` <-
function( Paz, Pdip,  Taz,  Tdip)
  {  
#### /**** convert P and T axis to cartesian coordinates  ***/
    DEG2RAD = pi/180
    RAD2DEG = 180/pi
    
    zp= cos(DEG2RAD*Pdip  );
    hp= sin(DEG2RAD*Pdip);
    yp= hp*cos(DEG2RAD*Paz );
    xp= hp*sin(DEG2RAD*Paz);
    
    P1 = list(x=xp, y = yp,z = zp) ;

    zt= cos(DEG2RAD*Tdip);
    ht= sin(DEG2RAD*Tdip);
    yt= ht*cos(DEG2RAD*Taz);
    xt= ht*sin(DEG2RAD*Taz);
    
    T1 = list(x=xt, y =yt, z = zt);
    
####/*** now get the cross product B =   P X T  ***/
    
    B = cross.prod(P1, T1);

    if(B$z>=1) B$z=1
    if(B$z<= -1 ) B$z= -1
    
    Bdip =RAD2DEG* acos(B$z);
    Baz= RAD2DEG*atan2( B$x, B$y);


    if(Baz<0.0) Baz=360+Baz;
    
    
    if(B$z>0.0)
      {
        Baz =    Baz-180;
        Bdip  =      90-Bdip;
      }
    else
      {
        Baz =    Baz;
        Bdip  =  90-Bdip;
      }
    
    if(Bdip<0) Bdip = -1.0*Bdip;

    return(list(Bdip=Bdip, Baz=Baz))
  }


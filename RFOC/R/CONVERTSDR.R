`CONVERTSDR` <-
function(strike, dip, rake)
{
 ### input is strike dip and rake in degrees
 ##  note in calculations we use the down dip directin here and not strike!
  ## dip is measured from the horizontal NOT from the  NADIR
  DEG2RAD = pi/180;
  RAD2DEG = 180/pi;
    dipdir = strike + 90.;
     signforp = 1;

    #  /* Compute F plane dip direction and dip */
    phif = RPMG::fmod(dipdir, 360.);
    deltaf = dip;

    #  /* Compute U slip vector */
    if (rake > 90.0) {  # /* this corrects for rakes past 90.0 by reflection */
	tmprake = rake - 180.;
	signforp = -1.;
    } else if (rake < -90.0) {
	tmprake = rake + 180.;
	signforp = -1.;
    } else {
	tmprake = rake;
	signforp = 1.;
    }
    
    temp1 = tmprake * DEG2RAD;
    temp2 = dip * DEG2RAD;
    deltau = asin(-1*sin(temp1) * sin(temp2) );      # /* leave in radians */
    tandy =  tan(deltau) / tan(temp2)
    if( tandy<(-1) ) tandy= (-1)
    if( tandy>(1) ) tandy= (1)
    
    phiu = dipdir - 90. + asin( tandy) *RAD2DEG;
    deltau = deltau*RAD2DEG;
    phiu = RPMG::fmod(phiu+360., 360.);
     # /* reflect(&phiu, &deltau); */

     # /* Compute G plane dip direction and dip */
    deltag = 90 - deltau;
    if (deltag < 0.) {
	phig = phiu;
	deltag = -deltag;
    } else if (deltag > 90.) {
	phig = phiu;
	deltag = 180. - deltag;
    } else {
	phig = phiu + 180.;
    }


    #################   this is not working?????
    if(deltaf==90 & rake==90)
      {
        ############   special case where above breaks down
       phig =  180 - phif 
       ###print(paste("special case", phig) )
        
      }
    
    
    phig = RPMG::fmod(phig, 360.);

     # /* Compute V slip vector */
    deltav = 90 - deltaf;
    phiv = RPMG::fmod(phif + 180, 360.);
     # /* reflect(&phiv, &deltav); */

     # /* Find P axis, as 1/2 way between U and V */
    U = TOCART.DIP(phiu, deltau);
    V = TOCART.DIP(phiv, deltav);
     # /* reverse second vector if necessary */
    
    V$x = signforp*V$x  ; V$y = signforp*V$y; V$z = signforp*V$z;
    
     # /* add together to form intermediate vector */
     
                                        #  P = TOSPHERE.DIP( U$x+V$x, U$y+V$y, U$z+V$z);

    
    P =  to.spherical( U$x+V$x, U$y+V$y, U$z+V$z)
    P$az = RPMG::fmod(P$az+ 360. , 360);
     
       P = REFLECT(P);

    ##############   here I ran into some problems
     ############   if P$dip>90  
    
     # /* Finally, find T axis orthogonal to P axis - this is easy */
     
     x2 = -V$x;
     y2 = -V$y;
     z2 = -V$z;
    T = to.spherical( U$x+x2, U$y+y2, U$z+z2 );
    T$az = RPMG::fmod( T$az , 360.);
     
    T = REFLECT(T);
    
     # /* Finally reflect U and V vectors */
    U = REFLECT(U);
    V = REFLECT(V);
     F = TOCART.DIP(phif, deltaf);
     G = TOCART.DIP(phig, deltag);

   M =list( az1=0, d1=0,  az2=0, d2=0, uaz=0, ud=0, vaz=0, vd=0, paz=0, pd =0, taz=0, td=0)
 
   M$az1=F$az;  M$d1=F$dip;
   M$az2=G$az;  M$d2=G$dip ;
   M$uaz=U$az;  M$ud=U$dip;
   M$vaz=V$az;  M$vd=V$dip;
   M$paz=P$az;  M$pd=P$dip;
   M$taz=T$az;  M$td=T$dip ;

  mc = structure( list(strike=strike, dipdir=dipdir, dip=dip, rake=rake, F=F, G=G, U=U, V=V, P=P, T=T, M=M), class="MEC" )
 
    return(mc);
}


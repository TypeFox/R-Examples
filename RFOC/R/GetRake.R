`GetRake` <-
function( az1, dip1, az2, dip2,  dir)
{
# az1=345; dip1=25; az2=122; dip2=71 ;

#   float *dipaz1,float *rake1,float *dipaz2,float *rake2
# /* 
# c J.C. Pechmann, July 1986
#     converted to C and modified  J.M.Lees Sept. 1994
# c
# c Rakcal calculates dip azimuths (dipaz1,dipaz2) and rake angles (rake1,rake2)
# c for a focal mechanism from the strike azimuths (az1,az2) and dips (dip1,dip2)
# c of the planes.  The strike azimuths must be in degrees measured clockwise from
# c north, with the nodal plane dipping down to the right of the strike direction.
# c Dip angles must be in degrees measured downward from the horizontal. Set 
# c dir= +1.0 if the faulting has a reverse component to it and set dir= -1.0
# c if the faulting has a normal component to it.  For strike-slip faulting, set
# c dir= +1.0 if the plane with the smaller strike azimuth is right lateral.
# c All angles returned by this subroutine are in degrees.
# c
# c
# c   
#   Calculate dip azimuths in degrees
# */

#  double rdip1,raz2, rdip2, rdipd1,rdipd2;
#  double zslip1, hslip1, yslip1, xslip1, zslip2;
#   double hslip2, yslip2, xslip2, zstrk1, ystrk1;
#   double raz1, xstrk1, zstrk2, ystrk2, xstrk2;
#    double dot1, dot2;     

DEG2RAD = pi/180
# print(paste(sep=" ", "in GetRake: plane 1", az1, dip1, "plane 2: ", az2, dip2,  "dir=",dir))
      dipaz1= az1 + 90.0;
#dipaz1= az1 - 90.0;
      if (dipaz1 >= 360.0) { dipaz1= RPMG::fmod(dipaz1, 360)  }
      dipaz2 = az2 + 90.0;
#dipaz2 = az2 - 90.0;
      if (dipaz2 >= 360.0) { dipaz2=  RPMG::fmod(dipaz2, 360.0) }

# /*     Convert angles to radians*/
      raz1=   az1*DEG2RAD;
      rdip1=  dip1*DEG2RAD;
      raz2=   az2*DEG2RAD;
      rdip2=  dip2*DEG2RAD;
      rdipd1= dipaz1*DEG2RAD;
      rdipd2= dipaz2*DEG2RAD;
# /*     Determine Cartersian coordinates for slip vectors (upper hemisphere)*/
      zslip1= cos(rdip2);
      hslip1= sin(rdip2);
      yslip1= hslip1*cos(rdipd2);
      xslip1= hslip1*sin(rdipd2);
      zslip2= cos(rdip1);
      hslip2= sin(rdip1);
      yslip2= hslip2*cos(rdipd1);
      xslip2= hslip2*sin(rdipd1);
#   
# /*     Determine Cartesian coordinates for unit vectors in the strike direction
# c     of each plane*/
      zstrk1= 0.0;
      ystrk1= cos(raz1);
      xstrk1= sin(raz1);
      zstrk2= 0.0;
      ystrk2= cos(raz2);
      xstrk2= sin(raz2);

# /*     Determine rake angles by taking dot products between slip vectors and
# c     strike vectors and then taking the inverse cosine*/

      dot1= xslip1*xstrk1 + yslip1*ystrk1 + zslip1*zstrk1;
      dot2= xslip2*xstrk2 + yslip2*ystrk2 + zslip2*zstrk2;

        if(dot1>1.0) { dot1 = 1.0; }

        if(dot1<(-1.0)) { dot1=-1.0; }

        if(dot2>1.0)  { dot2 = 1.0;}
        if(dot2<(-1.0) ) { dot2=-1.0;}

      rake1= acos(dot1)/DEG2RAD;
      rake2= acos(dot2)/DEG2RAD;
# /*     Adjust rake angles to match the Aki and Richards (p.106) sign convention.
# c     According to this convention, the rake angle is the angle between the
# c     direction of movement of the hanging wall and the strike direction of the
# c     footwall.  The sign is determined by the direction of movement of the
# c     hanging wall relative to the footwall, where up is defined as being 
# c     positive.  Thus, if the faulting has a reverse component to it, the rake
# c     angle is between 0 and 180 degrees.  If the faulting has a normal compo-
# c     nent to it, the angle is between 0 and -180 degrees.  For strike-slip 
# c     faulting, the hanging wall is defined as the right-hand block as viewed by
# c     an observer looking along the strike.  Thus, a left-lateral strike slip 
# c     fault has a rake angle of 0 and a right-lateral strike-slip fault has a
# c     rake angle of 180.*/
      if (dir<0) rake1= rake1- 180.0;
      if (dir<0) rake2= rake2- 180.0;
      if (rake1==(-180.0)) rake1= 180.0;
      if (rake2==(-180.0)) rake2= 180.0;






return(list(az1=az1, dip1=dip1, az2=az2, dip2=dip2,  dir=dir, rake1=rake1, dipaz1=dipaz1, rake2=rake2, dipaz2=dipaz2))
}


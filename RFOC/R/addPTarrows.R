`addPTarrows` <-
function(MEC)
{

   if(MEC$UP==TRUE)
    {
      az1=RPMG::fmod(MEC$P$az+180, 360)
    }
  else
    {
      az1=RPMG::fmod(MEC$P$az, 360)
    }
 
  dip = 90-MEC$P$dip
  Ppnt = focpoint(MEC$P$az, MEC$P$dip, pch=18, lab="P", UP=MEC$UP, PLOT=FALSE)
  JIMP  = RSEIS::TOCART(az1, dip)
   if(MEC$UP==FALSE ) JIMP$z = -JIMP$z
###  here x is pointing to the north, y to the east
#######   points( JIMP$y,JIMP$x) 
  len = 1
  headlen = len*.25
  basethick = len*0.02
  headlip = len*.008
  aglyph = Z3Darrow(len = len , basethick =basethick , headlen =headlen , headlip=headlip )
  Rview  =    ROTZ(0)%*%  ROTX(0)
### krappa = sqrt(2)
  krappa = sqrt(2)/2


  Vec1 = krappa*c(JIMP$y+Ppnt$x,JIMP$x+Ppnt$y,JIMP$z)
  Vec2 = c(0+Ppnt$x,0+Ppnt$y,0)
pcol = rgb(1,.5,.5)
  
  BOXarrows3D(Vec1[1], Vec1[2], Vec1[3],    Vec2[1], Vec2[2], Vec2[3] ,  aglyph=aglyph,  Rview=Rview, col=pcol)

  Ppnt = focpoint(MEC$P$az, MEC$P$dip, pch=18, lab="P", UP=MEC$UP, PLOT=TRUE)
#################################################
  Ppnt = focpoint(MEC$T$az, MEC$T$dip, pch=18, lab="P", UP=MEC$UP, PLOT=FALSE)
  if(MEC$UP==TRUE)
    {
      az1=RPMG::fmod(MEC$T$az+180, 360)
    }
  else
    {
      az1=RPMG::fmod(MEC$T$az, 360)
    }
  dip = 90-MEC$T$dip
  JIMP  = RSEIS::TOCART(az1, dip)
     if(MEC$UP==FALSE ) JIMP$z = -JIMP$z
###  here x is pointing to the north, y to the east
  len = 1
  headlen = len*.25
  basethick = len*0.02
  headlip = len*.008
  
  aglyph = Z3Darrow(len = len , basethick =basethick , headlen =headlen , headlip=headlip )

  Rview  =    ROTZ(0)%*%  ROTX(0)

  Vec2 = krappa*c(JIMP$y+Ppnt$x,JIMP$x+Ppnt$y,JIMP$z)
  Vec1 = c(0+Ppnt$x,0+Ppnt$y,0)

  tcol = rgb(.5,.5, 1)
  BOXarrows3D(Vec1[1], Vec1[2], Vec1[3],    Vec2[1], Vec2[2], Vec2[3] ,  aglyph=aglyph,  Rview=Rview, col=tcol)
  Ppnt = focpoint(MEC$T$az, MEC$T$dip, pch=18, lab="T", UP=MEC$UP, PLOT=TRUE)
}


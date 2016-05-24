`qpoint` <-
function(az, iang, col=2, pch=5, lab="", POS=4, UP=FALSE, PLOT=FALSE, cex=1 )
{
  #  plot a point on focal sphere,
  #  az = angle from north (degrees
  #  iang = angle from Z-down (not from horizontal)
  if(missing(POS)) { POS=4 ; }
  if(missing(PLOT)) { PLOT=TRUE }
  if(missing(cex)) { cex=1 }
   
  
  DEG2RAD = pi/180;
  sph = cos(DEG2RAD*iang);
  
  sph[sph>=0] = 0;
  sph[sph<0] =  1;
  
 
  if(missing(col))  {
    col=rep(1,length(az))
    col= col+sph
  }
  
  if(missing(pch))  {
    pch=rep(5,length(az))
    pch=pch+sph
  }
  
  if(missing(UP)) { UP=FALSE }


  if(UP==TRUE)
    {
      #  flip the orientation for upper sphere
      az = az+180
      # 
    }
  
  A = list(az=az, dip=iang)
   
  A$dip = iang
  A$az =az
  B = FixDip(A)

  
  trot =DEG2RAD* B$az
  tdip = B$dip

  xi =  DEG2RAD*tdip;

  tq = sqrt(2) * sin(xi / 2.0);
  
  pltx = tq * sin(trot);
  plty = tq * cos(trot);
  
  
  if(PLOT==TRUE)
    {
      points( pltx, plty , pch=pch, col=col, cex=cex)
      
      
      if(!missing(lab)) 
        {
          text( pltx, plty, labels=lab, pos=POS)
        }
    }
  return(list(x=pltx, y=plty))
}


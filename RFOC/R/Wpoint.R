`Wpoint` <-
function(az1, dip1, col=2, pch=5, lab="", UP=FALSE)
  {
  if(missing(col))  { col = 'red' }  
  if(missing(pch))  {  pch = 5 }  
  if(missing(lab))  {  lab=""  }
  if(missing(UP))  {UP=FALSE}


 DEG2RAD = pi/180;                                
  
  tdip = dip1
  if(UP==TRUE)
    {
      tdip=180 - dip1
    }
  else
    {
      tdip = 90-tdip
    }

  FLG = tdip>90

  tdip[FLG]=tdip[FLG]-90
  az1[FLG]=az1[FLG]-180

  if(missing(col))  {
    col=rep(2, length(az1))
    col[FLG]=3
   }
  if(missing(pch))  {
    pch=rep(5, length(az1))
    pch[FLG]=15
  }
  
  FLG2 = ( (dip1 == 0) | (dip1 == 180) )
          tdip[FLG2]=0
          az1[FLG2]=0
  
  trot =DEG2RAD*az1;
  xi =  DEG2RAD*tdip/2;
  tq =   tan(xi)
  pltx = tq * sin(trot);
  plty = tq * cos(trot);

  points( pltx, plty , pch=pch, col=col)


  if(!missing(lab)) 
    {
      text( pltx, plty, labels=lab, pos=4)
    }
  return(list(x=pltx, y=plty))


  }


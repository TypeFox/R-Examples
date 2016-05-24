`focpoint` <-
function(az1, dip1, col=2, pch=5, lab="", UP=FALSE, PLOT=TRUE )
{
                                        #    azimuth is degrees from north
                                        #    dip is degrees from down direction?  or degrees from horizontal?

  
                                        # print(paste(sep=" ", "focpoint", az1, dip1, lab)); 
  if(missing(col))  { col=rep(1,length(az1)) }
  if(missing(pch))  { pch=rep(5,length(az1)) }
  if(missing(PLOT)) { PLOT=TRUE; }


  if(length(col)<length(az1)) {  col=rep(col,length(az1))   }
  if(length(pch)<length(az1)) {  pch=rep(pch,length(az1))   }
  
                                        # if(missing(lab))  { lab="" }
  if(missing(UP)) { UP=FALSE }

  
  dip1 = 90-dip1 ;
  
  if(UP==TRUE)
    {
      az1=az1+180
      
    }
                                        # print(paste(sep=" ", "dip1=", dip1, "az1=", az1))
  DEG2RAD = pi/180;


  dflag = dip1>90
  
      dip1[dflag]=dip1[dflag]-90
      az1[dflag]=az1[dflag]-180
      pch[dflag]=15
      col[dflag]=3

  tdip = dip1;
  trot =DEG2RAD*az1;
  xi =  DEG2RAD*tdip;
  tq = sqrt(2) * sin(xi / 2.0);
  pltx = tq * sin(trot);
  plty = tq * cos(trot);
  if(PLOT==TRUE)
    {
      points( pltx, plty , pch=pch, col=col)
      
      if(!missing(lab)) 
        {
          text( pltx, plty, labels=lab, pos=4)
        }
    }
  invisible(list(x=pltx, y=plty))
}


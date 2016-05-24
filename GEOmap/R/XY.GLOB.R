`XY.GLOB` <-
function(x, y, PROJ.DATA)
  {

    if(PROJ.DATA$type==0)
      {
        
        LL = list(lon=x , lat=y)
        
      }
    if(PROJ.DATA$type==1)
      {
        LL = merc.sphr.ll(PROJ.DATA$LON0 , x , y)
      }
    if(PROJ.DATA$type==2)
      {
        LL = utm.sphr.ll( x , y, PROJ.DATA)
      }
    if(PROJ.DATA$type==3)
      {
        LL = lambert.cc.ll( x , y, PROJ.DATA)
      }
    if(PROJ.DATA$type==4)
      {
        LL = stereo.sphr.ll( x , y, PROJ.DATA)
      }
    if(PROJ.DATA$type==5)
      {
        LL = utm.elps.ll( x , y, PROJ.DATA)
      }
    if(PROJ.DATA$type==6)
      {
        LL = equid.cyl.ll(PROJ.DATA$LON0 ,PROJ.DATA$LAT0 , x , y)
      }
    if(PROJ.DATA$type==7)
      {
        LL = utm.wgs84.ll( x , y, PROJ.DATA)
      }
    
    if(PROJ.DATA$type==99)
      {
        LL = lcgc(PROJ.DATA$LAT0, PROJ.DATA$LON0,  x , y)
      }

    

    return(LL)

  }


`GLOB.XY` <-
function(LAT, LON, PROJ.DATA)
  {


    if(PROJ.DATA$type==0)
      {
        XY = list(x=LON, y=LAT)
      } 
    if(PROJ.DATA$type==1)
      {
        XY = merc.sphr.xy(PROJ.DATA$LON0 , LAT, LON)
      }
    if(PROJ.DATA$type==2)
      {
        XY = utm.sphr.xy( LAT, LON, PROJ.DATA)
      }
    if(PROJ.DATA$type==3)
      {
        XY = lambert.cc.xy( LAT, LON, PROJ.DATA)
      }
    if(PROJ.DATA$type==4)
      {
        XY = stereo.sphr.xy( LAT, LON, PROJ.DATA)
      }
    if(PROJ.DATA$type==5)
      {
        XY = utm.elps.xy( LAT, LON, PROJ.DATA)
      }
    if(PROJ.DATA$type==6)
      {
        XY = equid.cyl.xy(PROJ.DATA$LON0 ,PROJ.DATA$LAT0 , LAT, LON )
      }
    if(PROJ.DATA$type==7)
      {
        XY = utm.wgs84.xy( LAT, LON, PROJ.DATA)
      }



    
     if(PROJ.DATA$type==99)
      {
        XY =  gclc(PROJ.DATA$LAT0, PROJ.DATA$LON0, LAT, LON)
      }

    return(XY)

  }


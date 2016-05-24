KINOUT<-function(MAP, LLlim, projtype=2)
  {
    ###   return vecctor of strokes that are
    ###  in or out of a target region

    if(missing(projtype))  projtype=2
   
    if(!is.list(LLlim))
      {
### vlim = c( LLlim$lon[1],LLlim$lat[1], LLlim$lon[2],LLlim$lat[2])

        vlim = LLlim
        LLlim = list(lon=0, lat=0)
        LLlim$lon = vlim[c(1,3)]
        LLlim$lat =vlim[c(2,4)]
      }


    lons = RPMG::fmod(LLlim$lon, 360)
   ##   rlon = range(lons)

    if(lons[2]<lons[1])
      {
        lons[2] = lons[2]+360
      }

    lats = LLlim$lat
    cenlon = mean(lons)
    cenlat = mean(lats)

    
    PROJ = setPROJ(type=projtype, LAT0=cenlat, LON0=cenlon)

      XYLIM =  GLOB.XY(lats, lons , PROJ)
      LLlim$x = XYLIM$x
      LLlim$y = XYLIM$y


## print(XYLIM)


    
      y3 =XYLIM$y[1]
      y4 =XYLIM$y[2]
      x3 =XYLIM$x[1]
      x4 = XYLIM$x[2]
      
      XYMAP =  GLOB.XY( MAP$POINTS$lat, MAP$POINTS$lon, PROJ)
    
    IN = vector()
    kout = 0
    for(kin in 1:length(MAP$STROKES$num))
      {
        j1 = MAP$STROKES$index[kin]+1
        j2 = j1+MAP$STROKES$num[kin]-1
        
        if(j1>0 & j2>0 & j2-j1 >=0)
          {
            JEC = j1:j2
            x = XYMAP$x[JEC]
            y = XYMAP$y[JEC]
            
            if(any(y>y3 & y<y4 & x>x3 & x<x4))
              {
                kout = kout+1
                IN[kout] = kin
                
              }
          }
        
        }
    
    invisible(IN)
    
  }

`LOCPOLIMAP` <-
function(P, MAP)
  {
   ###  if(missing(shiftlon)) { shiftlon=0 }
    KAPPA = rep(0, times=length(MAP$STROKES$num))
    for(i in 1:length(MAP$STROKES$num))
      {
        ##  print(MAP[[i]]$x)
        j1 = MAP$STROKES$index[i]+1
        j2 = j1+MAP$STROKES$num[i]-1
        
        LONS = RPMG::fmod(MAP$POINTS$lon[j1:j2], 360)
        
        LATS = MAP$POINTS$lat[j1:j2]

        DL = abs(diff(LONS))

        if(any(DL>180))
          {
           ### print(paste(sep=' ', "flippin ",i,  MAP$STROKES$nam[i]))  
            LONS[LONS>180] = LONS[LONS>180] -360
            if(P$lon>180) { P$lon = P$lon -360 }
          }
        
        POK = list(x=LONS , y=LATS)

        J = inpoly(P$lon, P$lat, POK)
        KAPPA[i] = J
      ###  print(paste(sep=' ', J, MAP$STROKES$nam[i]))
      }

    return(KAPPA)
  }


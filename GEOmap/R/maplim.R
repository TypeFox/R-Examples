maplim<-function(lat, lon, pct=.1)
  {
###  pct is a percentage to expand (or shrink)
    
    rlat = range(lat, na.rm=TRUE)
    rlon = range(lon, na.rm=TRUE)

    if(pct!=0)
      {
        glat = expandbound(rlat, pct)
        glon = expandbound(rlon, pct)
      }
    else
      {
        glat = rlat
        glon = rlon

      }

    LON = RPMG::fmod(glon, 360)

    G = list(lat=glat, lon=glon, LON=LON, lim=c(glon[1], glat[1], glon[2], glat[2]), LIM=c(LON[1], glat[1], LON[2], glat[2]))

    return(G)
    
  }

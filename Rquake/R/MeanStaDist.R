MeanStaDist<-function(Ldat)
  {
    ########    calculate the mean km distance of
    ########   a set of Lat-lon pairs

    ww = which(!is.na(Ldat$lat) & !is.na(Ldat$lon) ) 
    lat = Ldat$lat[ww]
    lon = Ldat$lon[ww]
    if(length(lat)<1) return(NULL)
    
    MLAT = median(lat, na.rm = TRUE)
    MLON = median(lon)
    proj = GEOmap::setPROJ(type=2, LAT0=MLAT, LON0=MLON)
    XY = GEOmap::GLOB.XY(lat, lon, proj)
    DXY = dist(cbind(XY$x, XY$y))
    Mdistwt = mean(DXY)
    return(Mdistwt)
  }

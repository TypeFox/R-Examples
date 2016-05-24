
getGAP<-function( EQ, Ldat, PLOT=FALSE)
{
  
  MLAT = median(EQ$lat)
  MLON = median(EQ$lon)

 
  wsta = which( !duplicated(Ldat$name) )
  ###    Ldat$name[wsta]

  
  proj = GEOmap::setPROJ(type=2, LAT0=MLAT, LON0=MLON)
  XY = GEOmap::GLOB.XY(Ldat$lat[wsta],  Ldat$lon[wsta] , proj)
  eqxy = GEOmap::GLOB.XY(EQ$lat, EQ$lon, proj)

  gang = atan2(eqxy$y-XY$y ,  eqxy$x-XY$x)
### subtract off the first one
  gang = gang-gang[1]
###  make sure nothing crosses 2pi
  gang = RPMG::fmod(gang, 2*pi)

  ord1 = order(gang)

  rang = c( gang[ord1] , 2*pi)

  
  wm = which.max( abs( diff(rang) ) )

  gips = diff(rang)*180/pi
  gap = gips[wm]

  if(PLOT)
    {
      ord2 = c(ord1, 1)
      i1 = ord2[wm]
      i2 = ord2[wm+1]

      plot(XY$x, XY$y, asp=1)
      segments(eqxy$x, eqxy$y, XY$x, XY$y, col='blue')
      text(XY$x[ord1], XY$y[ord1], labels=1:length(XY$x), pos=1)
      segments( XY$x[i1]  ,  XY$y[i1], XY$x[i2]  ,  XY$y[i2], col='red')
     ##     ang1 = 180*atan2( XY$y[i1],  XY$x[i1])/pi
    ##      ang2 =  180*atan2( XY$y[i2],  XY$x[i2])/pi
    ##      rad = sqrt(XY$x[i1]^2  + XY$y[i1]^2 )
    ##   hh  = darc(rad = rad, ang1 =ang1 , ang2 =ang2 , x1 = 0, y1 = 0, n = 1)
      ##    lines(hh, col='green')
      
    }

  return(gap)
}


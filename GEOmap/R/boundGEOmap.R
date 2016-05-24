`boundGEOmap` <-
function(MAP, NEGLON=FALSE, projtype=2)
{
  if(missing(NEGLON)) { NEGLON=FALSE }
  if(missing(projtype)) { projtype=2 }

  
  ###############  if NEGLON=FALSE, convert all negative lons to positive
  ###########  if true, allow neg lons to stay negative
  ## caluclate the LL bounds of strokes in a GEOmap
  
 Kstroke = length(MAP$STROKES$num)

 LAT1 = rep(NA, length=Kstroke)
 LAT2 = rep(NA, length=Kstroke)
 LON1 = rep(NA, length=Kstroke)
 LON2 = rep(NA, length=Kstroke)
 
  for(i in 1:Kstroke)
    {
      j1 = MAP$STROKES$index[i]+1
      j2 = j1+MAP$STROKES$num[i]-1

      lon=MAP$POINTS$lon[j1:j2]
      lat=MAP$POINTS$lat[j1:j2]

      
      proj = setPROJ(type = projtype, LAT0 =lat[1] , LON0 =lon[1] )

      xy = GLOB.XY(lat, lon, proj)
      xr = range(xy$x)
      yr = range(xy$y)
      LL =  XY.GLOB(xr[1], yr[1], proj)
      UR =  XY.GLOB(xr[2], yr[2], proj)
      
      LAT1[i] = LL$lat
      LAT2[i] = UR$lat
      LON1[i] = LL$lon
      LON2[i] = UR$lon
      
    }

 MAP$STROKES$LAT1=LAT1
 MAP$STROKES$LAT2=LAT2
  
 MAP$STROKES$LON1=RPMG::fmod( LON1, 360)
 MAP$STROKES$LON2=RPMG::fmod(LON2, 360)

  if(NEGLON==FALSE)
    {

      MAP$STROKES$LON1[MAP$STROKES$LON1<0] = RPMG::fmod(MAP$STROKES$LON1[MAP$STROKES$LON1<0], 360)
      MAP$STROKES$LON2[MAP$STROKES$LON2<0] = RPMG::fmod(MAP$STROKES$LON2[MAP$STROKES$LON2<0], 360)
      MAP$POINTS$lon[MAP$POINTS$lon<0] = RPMG::fmod(MAP$POINTS$lon[MAP$POINTS$lon<0] , 360)

    }

 return(MAP)
}


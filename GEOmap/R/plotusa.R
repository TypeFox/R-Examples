`plotusa` <- function(USAmap, LATS=c(22,49.62741), LONS=c(229.29389,296.41803), add=FALSE)
{
 
  if(missing(LATS)) { LATS=c(22,49.62741) }
  if(missing(LONS)) { LONS=c(229.29389,296.41803) }
  if(missing(add)) { add=FALSE }

  USALL=list()
  USALL$lat=LATS
  USALL$lon=LONS
  PROJ = setPROJ(type = 2, LAT0 =mean(USALL$lat), LON0 = mean(USALL$lon) )

  plotGEOmapXY(USAmap, LIM= c(USALL$lon[1], USALL$lat[1], USALL$lon[2], USALL$lat[2]    )  , 
               PROJ=PROJ, add=add, shiftlon=0, axes=FALSE, ann=FALSE)
  
  invisible(list(PROJ=PROJ, USALL=USALL, LIM= c(USALL$lon[1], USALL$lat[1], USALL$lon[2], USALL$lat[2]    )))
  
}

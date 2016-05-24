`PlotArrowsOnStaticMap` <-structure(function#plots arrows or segments on map 
###This function plots/overlays arrows or segments on a map.
(
  MyMap, ##<< map image returned from e.g. \code{GetMap()} 
  lat0, ##<< latitude valuesof points FROM which to draw.
  lon0, ##<< longitude values of points FROM which to draw.
  lat1=lat0, ##<< latitude valuesof points TO which to draw.
  lon1=lon0, ##<< longitude values of points TO which to draw.
  TrueProj = TRUE, ##<< set to FALSE if you are willing to accept some degree of inaccuracy in the mapping. In that case, the coordinates of the image are in lat/lon and the user can simply overly points/lines/axis without worrying about projections
  FUN = arrows, ##<<, plotting function to use for overlay; typical choices would be \link{arrows} and \link{segments} 
  add = FALSE, ##<< start a new plot or add to an existing
  verbose = 0, ##<< level of verbosity 
  ... ##<< further arguments to be passed to \code{FUN}
){
  ##seealso<< \link{PlotOnStaticMap} \link{arrows}
  if (TrueProj){
    Rcoords0 <- LatLon2XY.centered(MyMap,lat= lat0,lon= lon0);   
    Rcoords1 <- LatLon2XY.centered(MyMap,lat= lat1,lon= lon1);
  }  else { #no transformtion
    Rcoords0 <- list(newY= lat0,newX= lon0);   
    Rcoords1 <- list(newY= lat1,newX= lon1);
  }
  
  if (!add) tmp <- PlotOnStaticMap(MyMap, TrueProj = TrueProj, verbose=verbose);
  FUN(x0=Rcoords0$newX, y0=Rcoords0$newY, x1=Rcoords1$newX, y1=Rcoords1$newY, ...);
### return value of \code{FUN}
}, ex = function(){
    	MyMap <- GetMap(center=c(lat=40.7,lon=-74), zoom=11)
	PlotArrowsOnStaticMap(MyMap, lat0=40.69, lon0=-73.9, lat1=40.71, lon1=-74.1, col = 'red')   

})




`TextOnStaticMap` <-structure(function#plots text on map
### TextOnStaticMap draws the strings given in the vector labels at the coordinates given by x and y on a map. y may be missing since xy.coords(x,y) is used for construction of the coordinates.
(
  MyMap, ##<<  map image returned from e.g. \code{GetMap()}
  lat, ##<< latitude where to put text.
  lon, ##<< longitude where to put text.
  labels = seq_along(lat), ##<< a character vector or \link{expression} specifying the text to be written. An attempt is made to coerce other language objects (names and calls) to expressions, and vectors and other classed objects to character vectors by \link{as.character}. If labels is longer than x and y, the coordinates are recycled to the length of labels.
  TrueProj = TRUE, ##<< set to FALSE if you are willing to accept some degree of inaccuracy in the mapping. In that case, the coordinates of the image are in lat/lon and the user can simply overly points/lines/axis without worrying about projections
  FUN = text, ##<< overlay function, typical choice would be \link{text}
  add = FALSE, ##<< start a new plot or add to an existing
  verbose = 0, ##<< level of verbosity 
  ... ##<< further arguments to be passed to \code{FUN}
){
      
  if (TrueProj){
    Rcoords <- LatLon2XY.centered(MyMap,lat= lat,lon= lon);   
  }  else { #no transformtion
    Rcoords <- list(newY= lat,newX= lon);   
  }
  
  if (!add) tmp <- PlotOnStaticMap(MyMap, TrueProj = TrueProj, verbose=verbose);
  FUN(x=Rcoords$newX, y=Rcoords$newY, labels = labels, ...);
### return value of \code{FUN}
}, ex = function(){
 lat = c(40.702147,40.718217,40.711614);
  lon = c(-74.012318,-74.015794,-73.998284);
  center = c(mean(lat), mean(lon));
  zoom <- min(MaxZoom(range(lat), range(lon)));
  
 
 MyMap <- GetMap(center=center, zoom=zoom,markers = paste0("&markers=color:blue|label:S|",
          "40.702147,-74.015794&markers=color:green|label:G|40.711614,-74.012318&markers=",
           "color:red|color:red|label:C|40.718217,-73.998284"), destfile = "MyTile1.png");
  TextOnStaticMap(MyMap, lat=40.711614,lon=-74.012318, "Some Text", cex=2, col = 'red')
  

})


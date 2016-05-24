`Tile2R` <-structure(function#simple utility to offset and scale XY coordinates with respect to the center
### simple utility to offset and scale XY coordinates with respect to the center
(
  points, ##<< XY coordinates returned by e.g. \link{LatLon2XY}
  center  ##<<  XY coordinates of center returned by e.g. \link{LatLon2XY}
){
##details<< mainly used for shrinking the size of a tile to the minimum size.

  X <- 256* (points$Tile[,"X"] - center$Tile[,"X"]) + (points$Coords[,"x"] - center$Coords[,"x"]) ;
  Y <- -256*(points$Tile[,"Y"] - center$Tile[,"Y"]) - (points$Coords[,"y"] - center$Coords[,"y"]) ;
  return(list(X=X,Y=Y))
### list with X and Y pixel values
}, ex = function(){
 latR <- c(34.5,34.9);
  lonR <- c(-100.3, -100);
  lat.center <- 34.7;
  lon.center <- -100.2;
  zoom = 10;
  ll <- LatLon2XY(latR[1], lonR[1], zoom);#lower left corner
  ur <- LatLon2XY(latR[2], lonR[2], zoom );#upper right corner
  cr <- LatLon2XY(lat.center, lon.center, zoom );#center
  ll.Rcoords <- Tile2R(ll, cr);
  ur.Rcoords <- Tile2R(ur, cr);

})



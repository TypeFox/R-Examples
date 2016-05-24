`XY2LatLon` <-structure(function#computes the centered coordinate transformation from lat/lon to map tile coordinates 
###The function XY2LatLon(MyMap, X,Y,zoom) computes the coordinate transformation from map tile coordinates to lat/lon given a map object.
(
  MyMap, ##<< map object
  X, ##<< latitude values to transform
  Y, ##<< longitude values to transform
  zoom ##<< optional zoom level. If missing, taken from \code{MyMap}
){
#X and Y are the centered integer pixel values, i.e. they should be zero in the center of the matrix/image !

   lat.center <- MyMap[[1]];
   lon.center <- MyMap[[2]];
   if (missing(zoom)) zoom <- MyMap[[3]];
      
   mycenter <- LatLon2XY(lat.center,lon.center,zoom);
   
   #first transform to original x,y coordinates
   x <- mycenter$Tile[,"X"] + (X+mycenter$Coords[,"x"])/256;
   y <- mycenter$Tile[,"Y"] - (Y-mycenter$Coords[,"y"])/256;
   
   ytilde <- 1 - y/2^(zoom-1);
   yy = (exp(2*pi* ytilde) - 1)/(exp(2*pi* ytilde) + 1);
   ShiftLat <- function(yy){
   	 n = c(-1,0,1);
     #print(180*asin( yy )/pi)
     #lat = 2*pi*(n+1) - asin( yy ) ; 
     lat = 2*pi*(n) + asin( yy ) ; 
     lat <- lat[which(lat <= pi/2 & lat > -pi/2 )];
     #convert back to degrees: 
     lat <- 180 * lat/ pi;
     return(lat);
   }
   lat <- sapply(yy, ShiftLat);
   #longitude is easy:
   lon =  180*( x/2^(zoom-1) - 1 );
##seealso<< \link{LatLon2XY} \link{Tile2R}

   return(cbind(lat=lat, lon=lon));
### properly scaled and centered (with respect to the center of \code{MyMap} ) coordinates  
###  \item{lon }{longitude}
###  \item{lat }{latitude}
}, ex = function(){
#quick test:

  zoom=12;MyMap <- list(40,-120,zoom, url="google");
  LatLon <- c(lat = 40.0123, lon = -120.0123);
  Rcoords <- LatLon2XY.centered(MyMap,LatLon["lat"],LatLon["lon"])
  newLatLon <- XY2LatLon(MyMap, Rcoords$newX, Rcoords$newY)
  max(abs(newLatLon - LatLon));

#more systematic:
 for (zoom in 2:10){
   cat("zoom: ", zoom, "\n");
   MyMap <- list(40,-120,zoom, url="google");
   LatLon <- c(lat = runif(1,-80,80), lon = runif(1,-170,170));
   Rcoords <- LatLon2XY.centered(MyMap,LatLon["lat"],LatLon["lon"])
   newLatLon <- XY2LatLon(MyMap, Rcoords$newX, Rcoords$newY)
   if(max(abs(newLatLon - LatLon)) > 0.0001) print(rbind(LatLon, newLatLon));
 }

})



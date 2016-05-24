`GetMap.bbox` <- structure(function
### Wrapper function for \link{GetMap}. Query the Google server for a static map tile, defined primarily by its lat/lon range and/or center and/or zoom. 
### Multiple additional arguments allow the user to customize the map tile.
( lonR, ##<< longitude range
  latR, ##<<  latitude range
  center, ##<< optional center
  size = c(640,640), ##<< desired size of the map tile image. defaults to maximum size returned by the Gogle server, which is 640x640 pixels 
  destfile = "MyTile.png", ##<<  File to load the map image from or save to, depending on \code{NEWMAP}.
  MINIMUMSIZE = FALSE, ##<< reduce the size of the map to its minimum size that still fits the lat/lon ranges ?
  RETURNIMAGE = TRUE, ##<< return image yes/no default: TRUE
  GRAYSCALE =FALSE, ##<< Boolean toggle; if TRUE the colored map tile is rendered into a black & white image, see \link{RGB2GRAY}
  NEWMAP = TRUE, ##<< if TRUE, query the Google server and save to \code{destfile}, if FALSE load from destfile.
  zoom,  ##<< Google maps zoom level. optional
  verbose=0, ##<< level of verbosity
  SCALE = 1, ##<< use the API's scale parameter to return higher-resolution map images. The scale value is multiplied with the size to determine the actual output size of the image in pixels, without changing the coverage area of the map
  ... ##<< extra arguments to \link{GetMap} 
){
  	if (missing(zoom)) zoom <- min(MaxZoom(latR, lonR, size));
 	if (missing(center)){
 	  lat.center <- mean(latR);#latR[1] + diff(latR)/2;
      lon.center <- mean(lonR);#lonR[1] + diff(lonR)/2;
    } else {
    	lat.center <- center[1];
    	lon.center <- center[2];
    }
     
    if (MINIMUMSIZE){
    	ll <- LatLon2XY(latR[1], lonR[1], zoom);#lower left corner
    	ur <- LatLon2XY(latR[2], lonR[2], zoom );#upper right corner
    	cr <- LatLon2XY(lat.center, lon.center, zoom );#center
    	ll.Rcoords <- Tile2R(ll, cr);
    	ur.Rcoords <- Tile2R(ur, cr);
    	if (verbose>1){
    		cat("ll:"); print(ll);print(ll.Rcoords)
    		cat("ur:"); print(ur);print(ur.Rcoords);
    		cat("cr:"); print(cr);
    	}
    	#return(list(ll,ur));
    	size[1] <- 2*max(c(ceiling(abs(ll.Rcoords$X)), ceiling(abs(ur.Rcoords$X)) ))+1;
    	size[2] <- 2*max(c(ceiling(abs(ll.Rcoords$Y)), ceiling(abs(ur.Rcoords$Y)) ))+1;
    	#size[1] <- ceiling(abs(256* (ur$Tile[,"X"] - ll$Tile[,"X"]) + (ur$Coords[,"x"] - ll$Coords[,"x"])));
    	#size[2] <- ceiling(abs(256* (ur$Tile[,"Y"] - ll$Tile[,"Y"]) + (ur$Coords[,"y"] - ll$Coords[,"y"])));
    	if (verbose) cat("new size: ", size, "\n")
    }
 	#if (NEWMAP) 
 	return(GetMap(center = c(lat.center, lon.center), zoom = zoom, size=size, destfile = destfile, RETURNIMAGE = RETURNIMAGE, GRAYSCALE = GRAYSCALE, SCALE=SCALE, verbose = verbose, ...));
 ### map tile

 }, ex = function(){
 	mymarkers <- cbind.data.frame(lat = c(38.898648,38.889112, 38.880940), 
          lon = c(-77.037692, -77.050273, -77.03660), size =  c('tiny','tiny','tiny'), 
          col = c('blue', 'green', 'red'), char = c('','',''));

	##get the bounding box:
  	bb <- qbbox(lat = mymarkers[,"lat"], lon = mymarkers[,"lon"]);
  
	##download the map:
  	MyMap <- GetMap.bbox(bb$lonR, bb$latR, destfile = "DC.png", GRAYSCALE =TRUE,
                markers = mymarkers);
 	##The function qbbox() basically computes a bounding box for the given lat,lon 
   #points with a few additional options such as quantile boxes, additional buffers, etc.  
  	bb <- qbbox(c(40.702147,40.711614,40.718217),c(-74.015794,-74.012318,-73.998284), 
            TYPE = "all", margin = list(m=rep(5,4), TYPE = c("perc", "abs")[1]));
 	##download the map:           
	MyMap <- GetMap.bbox(bb$lonR, bb$latR,destfile = "MyTile3.png", maptype = "satellite") 


})


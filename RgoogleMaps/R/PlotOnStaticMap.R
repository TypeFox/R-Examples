`PlotOnStaticMap` <-structure(function# overlays plot on background image of map tile
###This function is the workhorse of the package RgoogleMaps. It overlays plot on background image of map tile
(
  MyMap, ##<< optional map object 
  lat, ##<< latitude values to be overlaid
  lon, ##<< longitude values to be overlaid
  destfile, ##<< File to load the map image from or save to, depending on whether \code{MyMap} was passed. 
  zoom = NULL, ##<< Google maps zoom level. optional if \code{MyMap} is passed, required if not.
  size, ##<< desired size of the map tile image. defaults to maximum size returned by the Gogle server, which is 640x640 pixels  
  GRAYSCALE =FALSE, ##<< Boolean toggle; if TRUE the colored map tile is rendered into a black & white image, see \link{RGB2GRAY} 
  add=FALSE, ##<< start a new plot or add to an existing
  FUN = points, ##<< plotting function to use for overlay; typical choices would be \link{points} and \link{lines} 
  mar=c(0,0,0,0), ##<< outer margin in plot; if you want to see axes, change the default
  NEWMAP = TRUE, ##<< load map from file or get it "new" from the static map server
  TrueProj = TRUE, ##<< set to FALSE if you are willing to accept some degree of inaccuracy in the mapping. In that case, the coordinates of the image are in lat/lon and the user can simply overly points/lines/axis without worrying about projections
  axes= FALSE, ##<< overlay axes ?
  atX = NULL, ##<< numeric; position of ticks on x-axis; if missing, \link{axTicks} is called for nice values; see \link{axis}
  atY = NULL, ##<< numeric; position of ticks on y-axis; if missing, \link{axTicks} is called for nice values; see \link{axis}
  verbose = 0, ##<< level of verbosity 
  ... ##<< further arguments to be passed to \code{FUN}
){
  
  if (!NEWMAP & missing(MyMap) & missing(lat) & missing(lon) ){
    MyMap <- ReadMapTile(destfile, TRUE);
  } else if (missing(MyMap)) MyMap <- MapBackground(lat=lat,lon=lon, destfile =destfile, zoom=zoom, size=size, GRAYSCALE =GRAYSCALE, mar=mar, NEWMAP = NEWMAP, verbose=verbose);
  
  if (missing(size)) size <- dim(MyMap[[4]])[2:1];
  if (verbose == -1) browser();
  
  lat.center <- MyMap[[1]];
  lon.center <- MyMap[[2]];
  zoom <- MyMap[[3]];
  if (TrueProj & MyMap$url == "OSM") {
    print("Caution: map type is OpenStreetMap. Until we find the correct projection algorithm, we treat lat/lon like planar coordinates and set TrueProj = FALSE.")
    TrueProj = FALSE;
  }
  
  if ( !("BBOX" %in% names(MyMap)) ) MyMap$BBOX <- list(ll = XY2LatLon(MyMap, -size[1]/2 + 0.5, -size[2]/2 + 0.5), ur = XY2LatLon(MyMap, size[1]/2 - 0.5, size[2]/2 - 0.5) );
  
  if (verbose > 0) print(str(MyMap));
  if (!add) {
    #require(rimage)
    if ( all(is.na(charmatch(c("png","pdf","postscript", "jpeg"), names(dev.cur())))) ) {
      #if ( is.na(charmatch("null", names(dev.cur())) ) )  dev.off();
      #dev.new(width=7*size[1]/640, height=7*size[2]/640)
    }
    par(mar=mar);#par(pin=c(9,9))
    
    if (class(MyMap[[4]])[1] == "matrix"){
      x= seq(MyMap$BBOX$ll[2], MyMap$BBOX$ur[2], length= size[1]);
      y= seq(MyMap$BBOX$ll[1], MyMap$BBOX$ur[1], length= size[2]);
      image(x=x,y=y, z=MyMap[[4]], col = attr(MyMap[[4]], "COL"), xlab = "", ylab = "", axes = FALSE)
      #image(MyMap[[4]], "ind", col = attr(MyMap[[4]], "COL"))
    } else if (class(MyMap[[4]])[1] == "SpatialGridDataFrame"){
      image(MyMap[[4]], red=1, green=2, blue=3, axes = FALSE);
    } else if (class(MyMap[[4]])[1] %in% c("nativeRaster", "raster", "array")){
      if (exists("rasterImage")) { # can plot only in R 2.11.0 and higher
        plot(0:1,0:1,type="n", axes=FALSE, xlab="", ylab="", asp = size[2]/size[1]);
        tmp2 <- par('usr');
        updateusr(tmp2[1:2], x2=c(0,1), tmp2[3:4], y2=c(0,1) );
        #if (require(grid)) grid.raster(MyMap[[4]], width=1, height=1, y=0, just="bottom") else 
        rasterImage(MyMap[[4]], 0,0,1,1);
      } else {
        #myplot.imagematrix(x=x, y=y, z=MyMap[[4]], axes = FALSE);
      }
    } 
    tmp2 <- par("usr");
    offset = 1;
    if (TrueProj){
      SCALE = 2*MyMap$SCALE
      updateusr(tmp2[1:2], x2=c(-size[1]+offset, size[1]-offset)/SCALE, tmp2[3:4], y2=c(-size[2]+offset, size[2]-offset)/SCALE );
    } 
    if (axes){
      degreeAxis(1, MyMap=MyMap, at=atX); degreeAxis(2, MyMap=MyMap, at=atY);
    }
    #browser();
  }
  
  
  if (!missing(lat) & !missing(lon)){
    #if (TrueProj){
    Rcoords <- LatLon2XY.centered(MyMap,lat,lon);
    newX <- Rcoords$newX;
    newY <- Rcoords$newY;
    #browser()
    if (verbose) {
      print(range(newX, na.rm=TRUE));
      print(range(newY, na.rm=TRUE));
      #print(list(newX,newY));
    }
    #     } else {
    #        Rcoords<-transfXY(MyMap, lon, lat)
    #        newX <- Rcoords$newX;
    #        newY <- Rcoords$newY;
    #        #newX <- lon;
    #        #newY <- lat;
    #     }    
    
    FUN(newX, newY, ...)
  }
  
  #invisible(list(newX,newY)); 
  invisible(MyMap)
  ### the map object is returned via \code{invisible(MyMap)}
}, ex = function(){
#The first step naturally will be to download a static map from the Google server. A simple example:

  lat = c(40.702147,40.718217,40.711614);
  lon = c(-74.012318,-74.015794,-73.998284);
  center = c(mean(lat), mean(lon));
  zoom <- min(MaxZoom(range(lat), range(lon)));
  #this overhead is taken care of implicitly by GetMap.bbox();              
  MyMap <- GetMap(center=center, zoom=zoom,markers = paste0("&markers=color:blue|label:S|",
           "40.702147,-74.015794&markers=color:green|label:G|40.711614,-74.012318&markers=",
           "color:red|color:red|label:C|40.718217,-73.998284"), destfile = "MyTile1.png");
                 
   tmp <- PlotOnStaticMap(MyMap, lat = c(40.702147,40.711614,40.718217), 
                          lon = c(-74.015794,-74.012318,-73.998284), 
                          destfile = "MyTile1.png", cex=1.5,pch=20,
                          col=c('red', 'blue', 'green'), add=FALSE);
   #and add lines:
   PlotOnStaticMap(MyMap, lat = c(40.702147,40.711614,40.718217), 
                   lon = c(-74.015794,-74.012318,-73.998284), 
                   lwd=1.5,col=c('red', 'blue', 'green'), FUN = lines, add=TRUE)
   	

})




ReadMapTile <- structure(function#Read a bitmap image stored in the PNG format
###Reads an image from a PNG file/content into a raster array.
(
  destfile,  ##<< png file to read 
  METADATA = TRUE,  ##<< read MetaInfo as well ?
  native=TRUE ##<< determines the image representation - if FALSE then the result is an array, if TRUE then the result is a native raster representation, see \link{readPNG} in package png.
){
  fileBase <- substring(destfile,1, nchar(destfile)-4);
  fileExt <-  substring(destfile,nchar(destfile)-2,nchar(destfile));
  
 
  myTile <- readPNG(destfile, native=native);
 #if (exists("rasterImage")) myTile <- as.raster(myTile)

  #attr(myTile, "type") <- "rgb"; 
  size <- dim(myTile)[2:1];
  MetaInfo=NULL;
  if (METADATA) {
    try({
 	#load(paste(fileBase,"rda",sep="."));
 	VarsLoaded <- load(paste(destfile,"rda",sep="."));
 	if (is.element("MetaInfo", VarsLoaded))
 	  MyMap <- list(lat.center= MetaInfo$lat.center, lon.center=MetaInfo$lon.center, zoom=MetaInfo$zoom, 
                      myTile=myTile, BBOX = MetaInfo$BBOX, url = MetaInfo$url, size=size, SCALE=MetaInfo$SCALE);
 	return(MyMap);
   });
 } 
 return(myTile);
### map or tile object
})

	

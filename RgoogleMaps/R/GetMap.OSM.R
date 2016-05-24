`GetMap.OSM` <- structure(function#Query the Open Street Map server for map tiles instead of Google Maps
### The querying parameters for Open Street Maps are somewhat different in this version. 
### Instead of a zoom, center and size, the user supplies a scale parameter and a lat/lon bounding box. 
### The scale determines the image size.
(
  lonR=c(-74.02132,-73.98622), ##<< longitude range
  latR= c(40.69983,40.72595), ##<< latitude range
  scale= 20000, ##<< Open Street map scale parameter. The larger this value, the smaller the resulting map tile in memory. There is a balance to be struck between the lat/lon bounding box and the scale parameter.
  destfile = "MyTile.png", ##<<  File to load the map image from or save to, depending on \code{NEWMAP}.
  format = 'png', ##<< (optional) defines the format of the resulting image.
  RETURNIMAGE = TRUE, ##<< return image yes/no default: TRUE
  GRAYSCALE =FALSE, ##<< Boolean toggle; if TRUE the colored map tile is rendered into a black & white image, see \link{RGB2GRAY}
  NEWMAP = TRUE, ##<< if TRUE, query the Google server and save to \code{destfile}, if FALSE load from destfile.
  verbose=1, ##<< level of verbosity,
  ... ##<< extra arguments to be used in future versions
){
   ##note<<The OSM maptile server is frequently too busy to accomodate every request, so patience is warranted.
 	 options(scipen = 12);#to suppress scientific notation on the scale parameter
 	 
 	OSMbbox <- paste(lonR[1],latR[1],lonR[2],latR[2], sep=",")
 	#http://tile.openstreetmap.org/cgi-bin/export?bbox=-4.54,47.49,4.37,52.45&scale=14000000&format=png
    OSMurl <- 'http://tile.openstreetmap.org/cgi-bin/export?';
 	url <- paste(OSMurl, "bbox=", OSMbbox, "&scale=", scale, "&format=", format, sep="")
 	#OR, use zoom level (e.g. z=12 ):
 	# http://tah.openstreetmap.org/MapOf/index.php?long=-74.02132&lat=40.69983&z=12&w=256&h=256&format=png

 	if (verbose) print(url);
 	    
 	if (NEWMAP) ret <- download.file(url, destfile, mode="wb", quiet = FALSE);
    #BBOX <- list(ll = c(lonR[1], latR[1]), ur = c(lonR[2], latR[2]) );
    #thanks to Julien Barnier who corrected this bug:
    BBOX <- list(ll = c(latR[1], lonR[1]), ur = c(latR[2], lonR[2]))
	MetaInfo <- list(lat.center = mean(latR), lon.center  = mean(lonR), zoom = NULL, url = "OSM", BBOX = BBOX, scale=scale, SCALE = 1);
	save(MetaInfo, file = paste(destfile,"rda",sep="."));
 	if (verbose == -1) browser();

 	if (RETURNIMAGE){
 	  myMap <- ReadMapTile(destfile);
 	  if (GRAYSCALE) {
     	     myMap$myTile <- RGB2GRAY(myMap$myTile);
        }
	  return(myMap);	    
 	} else {
	  invisible(url)
	}
### map structure or URL used to download the tile.
 }, ex = function(){
if (interactive()) {
 	CologneMap <- GetMap.OSM(lonR= c(6.89, 7.09), latR = c(50.87, 51), scale = 150000, 
                            destfile = "Cologne.png");
	PlotOnStaticMap(CologneMap, mar=rep(4,4), NEWMAP = FALSE, TrueProj = FALSE, axes= TRUE);
		
	PrincetonMap <- GetMap.OSM(lonR= c(-74.67102, -74.63943), latR = c(40.33804,40.3556), 
                             scale = 12500, destfile = "Princeton.png");
	png("PrincetonWithAxes.png", 1004, 732)
      PlotOnStaticMap(PrincetonMap, axes = TRUE, mar = rep(4,4));
    dev.off()
 }
})



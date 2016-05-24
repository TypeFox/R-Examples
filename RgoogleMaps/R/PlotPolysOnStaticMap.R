`PlotPolysOnStaticMap` <-structure(function#plots polygons on map
###This function plots/overlays polygons on a map. Typically, the polygons originate from a shapefile.
(
  MyMap, ##<< map image returned from e.g. \code{GetMap()}
  polys, ##<<  polygons to overlay; these can be either of class \link[PBSmapping]{PolySet} from the package PBSmapping
### or of class \link[sp]{SpatialPolygons} from the package sp
  col, ##<< (optional) vector of colors, one for each polygon
  border = NULL, ##<< the color to draw the border. The default, NULL, means to use \link{par}("fg"). Use border = NA to omit borders, see \link{polygon}
  lwd = .25, ##<< line width, see \link{par}
  verbose = 0, ##<< level of verbosity 
  add=TRUE, ##<< start a new plot or add to an existing 
  textInPolys = NULL, ##<< text to be displayed inside polygons.
  ... ##<< further arguments passed to \code{PlotOnStaticMap}
){
  ##seealso<< \link{PlotOnStaticMap} \link{mypolygon}
  ##details<< 
  #print(str(polys))
  stopifnot(class(polys)[1] == "SpatialPolygons" | class(polys)[1] == "PolySet" | class(polys)[1] == "data.frame" | class(polys)[1] == "matrix")
  if (class(polys)[1] == "SpatialPolygons")
    polys = SpatialToPBS(polys)$xy
  
  Rcoords <- LatLon2XY.centered(MyMap,lat= polys[,"Y"],lon= polys[,"X"]);
  polys.XY <- as.data.frame(polys);
  polys.XY[,"X"] <- Rcoords$newX;
  polys.XY[,"Y"] <- Rcoords$newY;
  
  if ( !( "PID" %in% colnames(polys.XY)) ) 
     polys.XY[,"PID"] <- 1;
  if ( !( "SID" %in% colnames(polys.XY)) ) 
     polys.XY[,"SID"] <- 1;
  polys.XY[,"PIDSID"] <- apply(polys.XY[,c("PID","SID")],1,paste,collapse=":")
  if (!add) tmp <- PlotOnStaticMap(MyMap, verbose=0, ...)
  if (verbose>1) browser()
  #browser()
  if (!is.null(textInPolys)) Centers = PBSmapping::calcCentroid(polys.XY)
if (requireNamespace("PBSmapping", quietly = TRUE) & all(c("PID","X","Y","POS") %in% colnames(polys.XY)) ) {
	attr(polys.XY, "projection") <- NULL;
	usr <- par('usr')
	PBSmapping::addPolys(polys.XY,col=col, border = border, lwd = lwd, xlim =usr[1:2], ylim = usr[3:4],  ...)
	if (!is.null(textInPolys)) {
	  text(Centers[,"X"],Centers[,"Y"],textInPolys,cex=0.75, col = "blue")
	}
} else {
  if (!missing(col)) {
  	polys.XY[,"col"] <- col;
  	PIDtable <- as.numeric(table(polys.XY[,"PID"]));
  	SIDtable <- as.numeric(table(polys.XY[,"PIDSID"]))
  	if (length(SIDtable)==length(col))  polys.XY[,"col"] <- rep(col, SIDtable);
  	if (length(PIDtable)==length(col))  polys.XY[,"col"] <- rep(col, PIDtable);
  }
 
  if ( !( "col" %in% colnames(polys.XY)) )
     polys.XY[,"col"] <- rgb(.1,.1,.1,.05);
    #polys.XY[polys.XY[,"PID"] == 292,"col"] <- rgb(1,0,0,.75)
  #tmp <- by(polys.XY[,c("X","Y","col")], polys.XY[,"PID"], mypolygon, lwd = lwd, border = border);
  #if (0){
  	pids = unique(polys.XY[,"PIDSID"])
    
  	for (i in pids){
  	  jj = 	polys.XY[,"PIDSID"] == i;
  	  xx= polys.XY[jj,];
  	  if ( ( "POS" %in% colnames(xx)) ) xx <- xx[order(xx[,"POS"]),]
  	  polygon( xx[, c("X","Y")], col=xx[,"col"]);
  	  if (!is.null(textInPolys)) {
  	    text(Centers[i,"X"],Centers[i,"Y"],textInPolys[i], col = "blue")
  	  }
  	  #readLines(n=1)
  	}
  #}
  }
}, ex = function(){
if (interactive()){
  #require(PBSmapping);
  shpFile <- paste(system.file(package = "RgoogleMaps"), "/shapes/bg11_d00.shp", sep = "")
  #shpFile <- system.file('bg11_d00.shp', package = "RgoogleMaps");
  
  shp=importShapefile(shpFile,projection="LL");
  bb <- qbbox(lat = shp[,"Y"], lon = shp[,"X"]);
  MyMap <- GetMap.bbox(bb$lonR, bb$latR, destfile = "DC.png");
  PlotPolysOnStaticMap(MyMap, shp, lwd=.5, col = rgb(0.25,0.25,0.25,0.025), add = F);
  
  #Try an open street map:

  mapOSM <- GetMap.OSM(lonR=bb$lonR, latR=bb$latR, scale = 150000, destfile = "DC.png");
  PlotPolysOnStaticMap(mapOSM, shp, lwd=.5, col = rgb(0.75,0.25,0.25,0.15), add = F);

  #North Carolina SIDS data set:
  shpFile <- system.file("shapes/sids.shp", package="maptools");
  shp=importShapefile(shpFile,projection="LL");
  bb <- qbbox(lat = shp[,"Y"], lon = shp[,"X"]);
  MyMap <- GetMap.bbox(bb$lonR, bb$latR, destfile = "SIDS.png");
  #compute regularized SID rate
  sid <- 100*attr(shp, "PolyData")$SID74/(attr(shp, "PolyData")$BIR74+500)
  b <- as.integer(cut(sid, quantile(sid, seq(0,1,length=8)) ));
  b[is.na(b)] <- 1;
  opal <- col2rgb(grey.colors(7), alpha=TRUE)/255; opal["alpha",] <- 0.2;
  shp[,"col"] <- rgb(0.1,0.1,0.1,0.2);
  for (i in 1:length(b)) 
    shp[shp[,"PID"] == i,"col"] <- rgb(opal[1,b[i]],opal[2,b[i]],opal[3,b[i]],opal[4,b[i]]);
  PlotPolysOnStaticMap(MyMap, shp, lwd=.5, col = shp[,"col"], add = F);
  
  #compare the accuracy of this plot to a Google Map overlay:
  library(maptools);
  qk <- SpatialPointsDataFrame(as.data.frame(shp[, c("X","Y")]), as.data.frame(shp[, c("X","Y")]))
  sp::proj4string(qk) <- CRS("+proj=longlat");
  tf <- "NC.counties";
  SGqk <- GE_SpatialGrid(qk)
  png(file=paste(tf, ".png", sep=""), width=SGqk$width, height=SGqk$height,
  bg="transparent")
  par(mar=c(0,0,0,0), xaxs="i", yaxs="i");par(mai = rep(0,4))
  PBSmapping::plotPolys(shp, plt=NULL)
  dev.off()
  maptools::kmlOverlay(SGqk, paste(tf, ".kml", sep=""), paste(tf, ".png", sep=""));
  #This kml file can now be inspected in Google Earth or Google Maps
  
  #or choose an aspect ratio that corresponds better to North Carolina's elongated shape:
  MyMap <- GetMap.bbox(bb$lonR, bb$latR, destfile = "SIDS.png", size = c(640, 320), zoom = 7);
  PlotPolysOnStaticMap(MyMap, shp, lwd=.5, col = shp[,"col"], add = F);
 }
})



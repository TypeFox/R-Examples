# Purpose        : Generic method to plot time-series data in Google Earth 
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Dev Status     : Alpha
# Note           : plots are in the description tag;

setMethod("plotKML", "RasterBrickTimeSeries", function(
  obj,
  folder.name = normalizeFilename(deparse(substitute(obj, env=parent.frame()))),
  file.name = paste(folder.name, ".kml", sep=""),
  pngwidth = 680,
  pngheight = 180,
  pngpointsize = 14,
  kmz = get("kmz", envir = plotKML.opts),
  open.kml = TRUE,
  ...
){

  ## sampling locations:
  if(!("data" %in% slotNames(obj@sampled))){
    labs <- paste(obj@sampled@data[,1])
  } else {
    labs <- paste(1:length(obj@sampled))
  }
  ## Begin end times:
  TimeSpan.begin <- obj@TimeSpan.begin
  TimeSpan.end <- obj@TimeSpan.end
  ## copy mean times:
  obj@rasters <- setZ(obj@rasters, paste(as.POSIXct(unclass(as.POSIXct(TimeSpan.begin))+(unclass(as.POSIXct(TimeSpan.end))-unclass(as.POSIXct(TimeSpan.begin)))/2, origin="1970-01-01")))
  dtime = unclass(as.POSIXct(TimeSpan.end)) - unclass(as.POSIXct(TimeSpan.begin))

  ## open KML for writing:  
  kml_open(folder.name = folder.name, file.name = file.name)
  
  ## add a description for the whole folder:
  kml.out <- get("kml.out", envir=plotKML.fileIO)
  description_txt <- sprintf('<description>%s</description>', obj@rasters@title)
  parseXMLAndAdd(description_txt, parent=kml.out[["Document"]])  
  assign('kml.out', kml.out, envir=plotKML.fileIO)
  
  ## extract values at point locations:
  ov <- extract(obj@rasters, obj@sampled)
  png_names <- paste(obj@variable, "_timeseries_", 1:nrow(ov), ".png", sep="")
  html.table <- paste('<img src="', png_names, '" height="', pngheight, '" width="', pngwidth, '" align ="middle" />', sep = '')
  kml_layer.SpatialPoints(obj = obj@sampled, points_names = labs, html.table = html.table)
  
  ## plot rasters:
  kml_layer(obj = obj@rasters, dtime=dtime, ...) 

  ## plot the time-series data:
  for(i in 1:nrow(ov)){
    png(filename=png_names[i], width=pngwidth, height=pngheight, bg="white", pointsize=pngpointsize)
    par(mar=c(4.5,4.5,.8,.8))
    plot(as.Date(as.POSIXct(getZ(obj@rasters))), ov[i,], type="l", xlab="Date", ylab=obj@variable, col="grey", lwd=2)
    points(as.Date(as.POSIXct(getZ(obj@rasters))), ov[i,], pch="+", cex=.6)
    dev.off()
  }
  
  ## close the file:
  kml_close(file.name = file.name)
  if (kmz == TRUE){
      kml_compress(file.name = file.name)
  }
  ## open KML file in the default browser:
  if(open.kml==TRUE){
    kml_View(file.name)
  } else {
    message(paste("Object written to:", file.name))
  }
  
})


# end of script;
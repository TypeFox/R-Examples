# Purpose        : Writes a time series of rasters to a KML (all with the same legend);
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Pierre Roudier (pierre.roudier@landcare.nz); Dylan Beaudette (debeaudette@ucdavis.edu); 
# Status         : pre-alpha
# Note           : this function is only suitable for writing time-series of data i.e. multiple realizations of the same variables; we assume that the time dimension is set via the @zvalue slot;

kml_layer.RasterBrick <- function(
  obj,
  plot.legend = TRUE,  
  dtime = "", 
  tz = "GMT",
  z.lim = c(min(minValue(obj), na.rm=TRUE), max(maxValue(obj), na.rm=TRUE)),
  colour_scale = get("colour_scale_numeric", envir = plotKML.opts),
  home_url = get("home_url", envir = plotKML.opts),
  metadata = NULL,
  html.table = NULL,
  altitudeMode = "clampToGround",
  balloon = FALSE,
  png.width, 
  png.height,
  min.png.width = 800,
  png.type = "cairo-png",
  ...
  ){
   
  if(!is.numeric(obj@data@values)){
    stop('Values of class "numeric" required.') 
  }
  
  # Get our invisible file connection from custom environment
  kml.out <- get("kml.out", envir=plotKML.fileIO)
  
  # Checking the projection is geo
  prj.check <- check_projection(obj, control = TRUE)

  # Trying to reproject data if the check was not successful
  if (!prj.check) { obj <- reproject(obj) }

  # optional elevations:
  altitude <- charmatch("altitude", names(call))
  if(!is.na(altitude)){
    altitude <- eval(call[["altitude"]], nlayers(obj))
  } else {
    altitude <- rep(.all_kml_aesthetics[["altitude"]], length.out = nlayers(obj))
  }
  altitudeMode <- kml_altitude_mode(altitude)

  # Format the time slot for writing to KML:
  if(!any(class(getZ(obj)) %in% "POSIXct")|!any(class(getZ(obj)) %in% "character")){
    if(any(getZ(obj)=="")|is.null(getZ(obj))){
      obj <- setZ(obj, format(as.POSIXct(rev(as.Date(Sys.time())-1:nlayers(obj))), "%Y-%m-%dT%H:%M:%SZ"))
    }
      DateTime <- getZ(obj)[1:nlayers(obj)]
    }
    else { 
      DateTime <- getZ(obj)[1:nlayers(obj)] 
   }
  
  if(all(dtime==0)) {  
    when <- as.POSIXct(DateTime)
  }
  else {
    dtime <- mean(diff(unclass(as.POSIXct(DateTime))))    # estimate the time support (if not indicated)
    TimeSpan.begin <- format(as.POSIXct(unclass(as.POSIXct(DateTime)) - dtime/2, origin="1970-01-01"), "%Y-%m-%dT%H:%M:%SZ", tz=tz)
    TimeSpan.end <- format(as.POSIXct(unclass(as.POSIXct(DateTime)) + dtime/2, origin="1970-01-01"), "%Y-%m-%dT%H:%M:%SZ", tz=tz)
  }

  # Parse ATTRIBUTE TABLE (for each placemark):
  if(balloon & ("layernames" %in% slotNames(obj))){
      html.table <- .df2htmltable(data.frame(layernames=names(obj), zvalue=getZ(obj), unit=obj@unit))
  }

  # plot the legend (PNG)
  if(plot.legend == TRUE){
    legend_name <- paste(normalizeFilename(deparse(substitute(obj, env=parent.frame()))), "legend.png", sep="_")      
    colour_scale_legend <- colorRampPalette(colour_scale)(50)
    kml_legend.bar(x = z.lim, legend.file = legend_name, legend.pal = colour_scale_legend) 
  }

  message("Writing to KML...")
  # Name of the object
  pl1 = newXMLNode("Folder", parent=kml.out[["Document"]])
  pl2 <- newXMLNode("name", paste(class(obj)), parent=pl1)
  
  # Insert metadata:
  if(!is.null(metadata)){
    md.txt <- kml_metadata(metadata, asText = TRUE)
    txt <- sprintf('<description><![CDATA[%s]]></description>', md.txt)
    parseXMLAndAdd(txt, parent=pl1)
  }

  # Creating the PNG files using standard z.lim's:
  raster_name <- paste(normalizeFilename(names(obj)), ".png", sep="")

  # Plotting the image
  for(j in 1:length(raster_name)){
    if(missing(png.width)|missing(png.height)){
      png.width = ncol(raster(obj, j)); png.height = nrow(raster(obj, j))
    }
    ## minimum size of the images
    if(png.width<min.png.width){
      png.height <- round(min.png.width*png.height/png.width)
      png.width <- min.png.width 
    }
    png(filename = raster_name[j], bg = "transparent", type=png.type, width = png.width, height = png.height)
    par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
    colour_scale_legend <- colorRampPalette(colour_scale)(50)
    raster::image(raster(obj, j), col = colour_scale_legend, zlim = z.lim, frame.plot = FALSE, main="", maxpixels=ncell(raster(obj, j)))
    dev.off()
  }


  # Ground overlays:
  # =============
  if(length(html.table)>0 & all(dtime==0)){
    txtr <- sprintf('<GroundOverlay><name>%s</name><description><![CDATA[%s]]></description><TimeStamp><when>%s</when></TimeStamp><altitude>%.0f</altitude><altitudeMode>%s</altitudeMode><Icon><href>%s</href></Icon><LatLonBox><north>%.5f</north><south>%.5f</south><east>%.5f</east><west>%.5f</west></LatLonBox></GroundOverlay>', names(obj), html.table, when, altitude, rep(altitudeMode, length(raster_name)), paste(raster_name), rep(bbox(extent(obj))[2, 2], length(raster_name)), rep(bbox(extent(obj))[2, 1], length(raster_name)), rep(bbox(extent(obj))[1, 2], length(raster_name)), rep(bbox(extent(obj))[1, 1], length(raster_name))) 
  }
  else {
  if(length(html.table)>0 & any(!dtime==0)){  # with attributes / block temporal support 
    txtr <- sprintf('<GroundOverlay><name>%s</name><description><![CDATA[%s]]></description><TimeSpan><begin>%s</begin><end>%s</end></TimeSpan><altitude>%.0f</altitude><altitudeMode>%s</altitudeMode><Icon><href>%s</href></Icon><LatLonBox><north>%.5f</north><south>%.5f</south><east>%.5f</east><west>%.5f</west></LatLonBox></GroundOverlay>', names(obj), html.table, TimeSpan.begin, TimeSpan.end, rep(altitude, length(raster_name)), altitude, paste(raster_name), rep(bbox(extent(obj))[2, 2], length(raster_name)), rep(bbox(extent(obj))[2, 1], length(raster_name)), rep(bbox(extent(obj))[1, 2], length(raster_name)), rep(bbox(extent(obj))[1, 1], length(raster_name)))
  }
  else {
  if(is.null(html.table) & any(!dtime==0)){   # no attributes / block temporal support 
    txtr <- sprintf('<GroundOverlay><name>%s</name><TimeSpan><begin>%s</begin><end>%s</end></TimeSpan><altitude>%.0f</altitude><altitudeMode>%s</altitudeMode><Icon><href>%s</href></Icon><LatLonBox><north>%.5f</north><south>%.5f</south><east>%.5f</east><west>%.5f</west></LatLonBox></GroundOverlay>', names(obj), TimeSpan.begin, TimeSpan.end, altitude, rep(altitudeMode, length(raster_name)), paste(raster_name), rep(bbox(extent(obj))[2, 2], length(raster_name)), rep(bbox(extent(obj))[2, 1], length(raster_name)), rep(bbox(extent(obj))[1, 2], length(raster_name)), rep(bbox(extent(obj))[1, 1], length(raster_name)))
  }
  else {  # no attributes / point temporal support 
     txtr <- sprintf('<GroundOverlay><name>%s</name><TimeStamp><when>%s</when></TimeStamp><altitude>%.0f</altitude><altitudeMode>%s</altitudeMode><Icon><href>%s</href></Icon><LatLonBox><north>%.5f</north><south>%.5f</south><east>%.5f</east><west>%.5f</west></LatLonBox></GroundOverlay>', names(obj), when, altitude, rep(altitudeMode, length(raster_name)), paste(raster_name), rep(bbox(extent(obj))[2, 2], length(raster_name)), rep(bbox(extent(obj))[2, 1], length(raster_name)), rep(bbox(extent(obj))[1, 2], length(raster_name)), rep(bbox(extent(obj))[1, 1], length(raster_name)))
  }}}

  parseXMLAndAdd(txtr, parent=pl1)

  # Legend
  # ======================
  if(plot.legend == TRUE){
    txtso <- sprintf('<ScreenOverlay><name>Legend</name><Icon><href>%s</href></Icon><overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/><screenXY x="0" y="1" xunits="fraction" yunits="fraction"/></ScreenOverlay>', legend_name)
    parseXMLAndAdd(txtso, parent=kml.out[["Document"]])
  }
  
  # save results: 
  assign("kml.out", kml.out, envir=plotKML.fileIO)
  
}

setMethod("kml_layer", "RasterBrick", kml_layer.RasterBrick)

# end of script;
# Purpose        : Write a SpatialPixels object to KML;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : ; 
# Status         : pre-alpha
# Note           : ;

kml_layer.SpatialPixels <- function(  
  obj,
  subfolder.name = paste(class(obj)),
  raster_name,
  plot.legend = TRUE,
  metadata = NULL,
  png.width = gridparameters(obj)[1,"cells.dim"], 
  png.height = gridparameters(obj)[2,"cells.dim"],
  min.png.width = 800,
  TimeSpan.begin,
  TimeSpan.end,
  layer.name,
  png.type = "cairo-png",
  ...
  ){

  ## get our invisible file connection from custom evnrionment
  kml.out <- get("kml.out", envir=plotKML.fileIO)

  ## Checking the projection 
  prj.check <- check_projection(obj, control = TRUE)

  ## Parsing the call:
  call.lst <- substitute(list(...))
  call.lst <- as.list(call.lst)[-1]

  ## Check if any attribute has been selected:
  if (is.na(charmatch("colour", names(call.lst)))){
    stop("No attribute selected. Please use the colour = ... option.")
  }
  if(missing(layer.name)){
    layer.name = deparse(call.lst[["colour"]])
  }

  if(is.call(call.lst[["colour"]])|is.name(call.lst[["colour"]])){
    x <- data.frame(eval(call.lst[["colour"]], obj@data))
    names(x) <- deparse(call.lst[["colour"]])
    obj@data <- x
  } else { 
  if(is.numeric(call.lst[["colour"]])) {
    i_layer <- call.lst[["colour"]]
    if (nlayers(obj) > 1) {
      obj <- obj[i_layer]
    }
  } else { 
  if(is.character(call.lst[["colour"]])) {
    i_layer <- which(names(obj) == call.lst[["colour"]])
    if (nlayers(obj) > 1) {
      obj <- obj[i_layer]
    }
  }}}

  ## Trying to reproject data if the check was not successful
  if(!prj.check) {  
    obj <- reproject(obj) 
  }

  ## TH: this needs to be fixed
  r <- raster(obj)
  altitude <- charmatch("altitude", names(call.lst))
  if(!is.na(altitude)){
    altitude <- eval(call.lst[["altitude"]], nlayers(r))
  } else {
    altitude <- rep(.all_kml_aesthetics[["altitude"]], length.out = nlayers(r))
  }
  altitudeMode <- kml_altitude_mode(altitude, GroundOverlay=TRUE) 

  ## prepare the palette:
  if (!is.na(charmatch("colour_scale", names(call.lst)))){
    pal <- eval(call.lst[["colour_scale"]])
  } else {
  ## default colour palettes
    if (!is.factor(obj@data[,1])){
      pal <- get("colour_scale_numeric", envir = plotKML.opts)
    } else {
      pal <- get("colour_scale_factor", envir = plotKML.opts)
    }
  }

  if(is.factor(obj@data[,1])){
    colour_scale <- colorRampPalette(pal)(length(levels(obj@data[,1])))
  } else {
    colour_scale <- colorRampPalette(pal)(100)  
  }

  ## Transparency
  alpha <- charmatch("alpha", names(call.lst))
  if (!is.na(alpha)) {
    colour_scale <- kml_alpha(obj, alpha = eval(call.lst[["alpha"]], obj@data), colours = colour_scale, RGBA = TRUE)
  }

  ## Creating the PNG file
  if(missing(raster_name)){
    raster_name <- paste(normalizeFilename(as.character(call.lst[["colour"]])), ".png", sep="")
  }

  ## Plotting the image
  if(png.width<min.png.width){
     png.height <- round(min.png.width*png.height/png.width)
     png.width <- min.png.width  
  }
  png(filename = raster_name, bg = "transparent", type = png.type, width = png.width, height = png.height)
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  if(!is.na(charmatch("z.lim", names(call.lst)))){ 
    z.lim <- eval(call.lst[["z.lim"]]) 
    r <- calc(r, fun=function(x){ x[x < z.lim[1]] <- z.lim[1]; return(x)}) 
    r <- calc(r, fun=function(x){ x[x > z.lim[2]] <- z.lim[2]; return(x)})
    raster::image(r, col = colour_scale, zlim = z.lim, frame.plot = FALSE, main="", maxpixels=ncell(r))
  } else {
    if(is.factor(obj@data[,1])){ 
      f.breaks <- seq(0.5, length(levels(obj@data[,1]))+0.5)
      raster::image(r, col = colour_scale, breaks=f.breaks, frame.plot = FALSE, main="", maxpixels=ncell(r))
    } else {
      raster::image(r, col = colour_scale, frame.plot = FALSE, main="", maxpixels=ncell(r))
    }
  }
  dev.off()

  ## There is a bug in Google Earth that does not allow transparency of PNGs:
  ## http://groups.google.com/group/earth-free/browse_thread/thread/1cd6bc29a2b6eb76/62724be63547fab7
  ## Solution: add transparency using ImageMagick:
  convert <- get("convert", envir = plotKML.opts)
  if(nchar(convert)==0){
    plotKML.env(silent = FALSE, show.env = FALSE)
    convert <- get("convert", envir = plotKML.opts)
  } else {
    ## if it does manages to find ImageMagick:
    if(!nchar(convert)==0){
      system(paste(convert, ' ', raster_name, ' -matte -transparent "#FFFFFF" ', raster_name, sep=""))
    } else {
    warning("PNG transparency possibly ineffective. Install ImageMagick and add to PATH. See ?kml_layer.SpatialPixels for more info.")
    }
  }

  ## plot the legend (PNG)
  if(plot.legend == TRUE){
    if(missing(raster_name)){
      legend_name <- paste(as.character(call.lst[["colour"]]), "legend.png", sep="_")
    } else {
      legend_name <- paste(strsplit(raster_name, "\\.")[[1]][1], "legend.png", sep="_")  
    }
    if(!is.na(charmatch("z.lim", names(call.lst)))){
      kml_legend.bar(x = obj@data[,1], legend.file = legend_name, legend.pal = colour_scale, z.lim = eval(call.lst[["z.lim"]]))
    } else {
      kml_legend.bar(x = obj@data[,1], legend.file = legend_name, legend.pal = colour_scale)
    }
  }

  message("Writing to KML...")
  ## Folder name
  pl1 = newXMLNode("Folder", parent=kml.out[["Document"]])
  pl2 <- newXMLNode("name", subfolder.name, parent = pl1)

  ## Add time stamp if available:
  if(!missing(TimeSpan.begin)&!missing(TimeSpan.end)){
    if(TimeSpan.end<TimeSpan.begin){
      stop("Missing or invalid 'TimeSpan.begin' and/or 'TimeSpan.end'. See also 'kml_layer.STIDF'")
    } 
    TimeSpan.begin = format(TimeSpan.begin, "%Y-%m-%dT%H:%M:%SZ")
    TimeSpan.end = format(TimeSpan.end, "%Y-%m-%dT%H:%M:%SZ")
  }

  ## Insert metadata:
  if(!is.null(metadata)){
    md.txt <- kml_metadata(metadata, asText = TRUE)
    txt <- sprintf('<description><![CDATA[%s]]></description>', md.txt)
    parseXMLAndAdd(txt, parent=pl1)
  }

  ## Ground overlay
  ## =====================
  pl2b <- newXMLNode("GroundOverlay", parent = pl1)
  ## Creating a SpatialPixelsDataFrame object to be plotted
  pl3 <- newXMLNode("name", layer.name, parent = pl2b)
  if(!missing(TimeSpan.begin)&!missing(TimeSpan.end)){
    pl3a <- newXMLNode("TimeSpan", parent = pl2b)
    pl3a1 <- newXMLNode("begin", TimeSpan.begin, parent = pl3a)
    pl3a2 <- newXMLNode("end", TimeSpan.end, parent = pl3a)
  }
  pl3b <- newXMLNode("altitude", signif(altitude, 4), parent = pl2b)
  pl3b <- newXMLNode("altitudeMode", altitudeMode, parent = pl2b)
  pl3c <- newXMLNode("Icon", parent = pl2b)
  pl4 <- newXMLNode("href", raster_name, parent = pl3c)
  pl3d <- newXMLNode("LatLonBox", parent = pl2b)
  pl4b <- newXMLNode("north", bbox(obj)[2, 2], parent = pl3d)
  pl4c <- newXMLNode("south", bbox(obj)[2, 1], parent = pl3d)
  pl4d <- newXMLNode("east", bbox(obj)[1, 2], parent = pl3d)
  pl4e <- newXMLNode("west", bbox(obj)[1, 1], parent = pl3d)
  
  ## Legend
  ## ======================
  if(plot.legend == TRUE){
    txtso <- sprintf('<ScreenOverlay><name>Legend</name><Icon><href>%s</href></Icon><overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/><screenXY x="0" y="1" xunits="fraction" yunits="fraction"/></ScreenOverlay>', legend_name)
    parseXMLAndAdd(txtso, parent=kml.out[["Document"]])
  }
  
  ## save results: 
  assign("kml.out", kml.out, envir=plotKML.fileIO)

}

setMethod("kml_layer", "SpatialPixels", kml_layer.SpatialPixels)
setMethod("kml_layer", "SpatialGrid", kml_layer.SpatialPixels)

# end of script;
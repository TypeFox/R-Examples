# Purpose        : Write a raster layer to KML;
# Maintainer     : Pierre Roudier (pierre.roudier@landcare.nz);
# Contributions  : Tomislav Hengl (tom.hengl@wur.nl); Dylan Beaudette (debeaudette@ucdavis.edu); 
# Status         : pre-alpha
# Note           : Rasters can also be written as polygons; see "?grid2poly";

kml_layer.Raster <- function(
  obj,
  subfolder.name = paste(class(obj)),  
  plot.legend = TRUE,
  metadata = NULL,
  raster_name,
  png.width = ncol(obj), 
  png.height = nrow(obj),
  min.png.width = 800,
  TimeSpan.begin,
  TimeSpan.end,
  layer.name,
  png.type = "cairo-png",
  ...
  ){

  # get our invisible file connection from custom environment
  kml.out <- get("kml.out", envir=plotKML.fileIO)

  # Checking the projection 
  prj.check <- check_projection(obj, control = TRUE)

  # Parsing the call:
  call.lst <- substitute(list(...))
  call.lst <- as.list(call.lst)[-1]

  # Check if any attribute has been selected:
  if (is.na(charmatch("colour", names(call.lst)))){
    stop("No attribute selected. Please use the colour = ... option.")
  }
  if(missing(layer.name)){
    layer.name = deparse(call.lst[["colour"]])
  }

  if(is.call(call.lst[["colour"]])|is.name(call.lst[["colour"]])){
    x <- data.frame(getValues(obj))
    names(x) <- names(obj)
    x <- eval(call.lst[["colour"]], x)
    obj <- raster(obj)
    values(obj) <- x
  } else { 
  if(is.numeric(call.lst[["colour"]])) {
    i_layer <- call.lst[["colour"]]
    if (nlayers(obj) > 1) {
      obj <- raster(obj, layer = i_layer)
    }
  } else { 
  if(is.character(call.lst[["colour"]])) {
    i_layer <- which(names(obj) == call.lst[["colour"]])
    if (nlayers(obj) > 1) {
      obj <- raster(obj, layer = i_layer)
    }
  }}}

  # TH: this needs to be fixed
  altitude <- charmatch("altitude", names(call.lst))
  if(!is.na(altitude)){
    altitude <- eval(call.lst[["altitude"]], nlayers(obj))
  } else {
    altitude <- rep(.all_kml_aesthetics[["altitude"]], length.out = nlayers(obj))
  }
  altitudeMode <- kml_altitude_mode(altitude, GroundOverlay=TRUE) 

  # prepare the palette:
  if (!is.na(charmatch("colour_scale", names(call.lst)))){
    pal <- eval(call.lst[["colour_scale"]])
  } else {
  ## default colour palettes
    if (!is.factor(obj)){
      pal <- get("colour_scale_numeric", envir = plotKML.opts)
    } else {
      pal <- get("colour_scale_factor", envir = plotKML.opts)
    }
  }

  # Trying to reproject data if the check was not successful
  if(!prj.check) {  
    obj <- reproject(obj) 
  }

  if(is.factor(obj)){
    if(length(labels(obj))==0){
      warning("RasterLayer of type factor missing labels")
      colour_scale <- colorRampPalette(pal)(length(levels(as.factor(getValues(obj)))))      
    } else {
      colour_scale <- colorRampPalette(pal)(length(labels(obj)[[1]]))
    }
  } else {
    colour_scale <- colorRampPalette(pal)(100)  
  }

  # Transparency
  alpha <- charmatch("alpha", names(call.lst))
  if (!is.na(alpha)) {
    ## - constant transparency
    ## - raster index if a Stack
    ## - name of a layer if a Stack
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
  png(filename = raster_name, bg = "transparent", type=png.type, width=png.width, height=png.height)
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  if(!is.na(charmatch("z.lim", names(call.lst)))){ 
    z.lim <- eval(call.lst[["z.lim"]])
    obj <- calc(obj, fun=function(x){ x[x < z.lim[1]] <- z.lim[1]; return(x)}) 
    obj <- calc(obj, fun=function(x){ x[x > z.lim[2]] <- z.lim[2]; return(x)})
    raster::image(obj, col = colour_scale, zlim = z.lim, frame.plot = FALSE, main="", maxpixels=ncell(obj)) 
  } else {
    if(is.factor(obj)&!length(levels(obj))==0){ 
        f.breaks <- seq(0.5, length(levels(obj))+0.5)
        raster::image(obj, col = colour_scale, breaks=f.breaks, frame.plot = FALSE, main="", maxpixels=ncell(obj))
      } else {
        raster::image(obj, col = colour_scale, frame.plot = FALSE, main="", maxpixels=ncell(obj))
      }
  }
  dev.off()

  ## There is a bug in Google Earth that does not allow transparency of PNGs:
  # http://groups.google.com/group/earth-free/browse_thread/thread/1cd6bc29a2b6eb76/62724be63547fab7
  # Solution: add transparency using ImageMagick:
  convert <- get("convert", envir = plotKML.opts)
  if(nchar(convert)==0){
    plotKML.env(silent = FALSE, show.env = FALSE)
    convert <- get("convert", envir = plotKML.opts)
  } else {
    # if it does manages to find ImageMagick:
    if(!nchar(convert)==0){
      system(paste(convert, ' ', raster_name, ' -matte -transparent "#FFFFFF" ', raster_name, sep=""))
    } else {
    warning("PNG transparency possibly ineffective. Install ImageMagick and add to PATH. See ?kml_layer.Raster for more info.")
    }
  }

  ## plot the legend (PNG)
  if(plot.legend == TRUE){
    if(missing(raster_name)){
      legend_name <- paste(as.character(call.lst[["colour"]]), "legend.png", sep="_")
    } else {
      legend_name <- paste(strsplit(raster_name, "\\.")[[1]][1], "legend.png", sep="_")      
    }
    if(is.factor(obj)){
      x <- as.factor(getValues(obj))
      if(length(labels(obj))==0){
        levels(x) <- levels(as.factor(getValues(obj)))
      } else {
        levels(x) = labels(obj)[[1]]
      }      
      kml_legend.bar(x = x, legend.file = legend_name, legend.pal = colour_scale)
    } else {
      if(!is.na(charmatch("z.lim", names(call.lst)))){
        kml_legend.bar(x = getValues(obj), legend.file = legend_name, legend.pal = colour_scale, z.lim = eval(call.lst[["z.lim"]]))
       } else {
        kml_legend.bar(x = getValues(obj), legend.file = legend_name, legend.pal = colour_scale)
       } 
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
  pl4b <- newXMLNode("north", bbox(extent(obj))[2, 2], parent = pl3d)
  pl4c <- newXMLNode("south", bbox(extent(obj))[2, 1], parent = pl3d)
  pl4d <- newXMLNode("east", bbox(extent(obj))[1, 2], parent = pl3d)
  pl4e <- newXMLNode("west", bbox(extent(obj))[1, 1], parent = pl3d)
  
  ## Legend
  ## ======================
  if(plot.legend == TRUE){
    txtso <- sprintf('<ScreenOverlay><name>Legend</name><Icon><href>%s</href></Icon><overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/><screenXY x="0" y="1" xunits="fraction" yunits="fraction"/></ScreenOverlay>', legend_name)
    parseXMLAndAdd(txtso, parent=kml.out[["Document"]])
  }
  
  ## save results: 
  assign("kml.out", kml.out, envir=plotKML.fileIO)
  
}

setMethod("kml_layer", "RasterLayer", kml_layer.Raster)
setMethod("kml_layer", "RasterStack", kml_layer.Raster)

# end of script;

# Purpose        : Generic method to plot simulated equiprobable vectors in Google Earth 
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Dev Status     : Alpha
# Note           : it basically requires only a single input object;


## plot object:
setMethod("plotKML", "SpatialVectorsSimulations", function(
  obj,
  folder.name = normalizeFilename(deparse(substitute(obj, env=parent.frame()))),
  file.name = paste(folder.name, ".kml", sep=""),
  colour,
  grid2poly = FALSE,
  obj.summary = TRUE,
  plot.svar = FALSE,
  kmz = get("kmz", envir = plotKML.opts),
  open.kml = TRUE,
  ...
){
  
  # if missing colour, pick the first var on the list
  if(missing(colour)){ 
    obj@summaries@data[,"colour"] <- obj@summaries@data[,1] 
    message("Plotting the first variable on the list")  
  } else {
    if(is.name(colour)|is.call(colour)){
      obj@summaries@data[,"colour"] <- eval(colour, obj@summaries@data)
    } else {
      obj@summaries@data[,"colour"] <- obj@summaries@data[,as.character(colour)]  
    }
  }
  
  # error map (this assumes that it is always the 2nd on the list):
  colour.sd <- obj@summaries@data[,2]
  # mask out 0 pixels
  obj@summaries@data[,"colour"] <- ifelse(obj@summaries@data[,"colour"]==0, NA, obj@summaries@data[,"colour"])
  N.r <- length(obj@realizations)

  # summary properties of the RK model:
  if(obj.summary==TRUE){
    sel <- obj@summaries@data[,"colour"]>0
    md <- data.frame(Names=c("N.realizations", "avg.probability", "N.pixels"), Values=c(N.r, signif(mean(obj@summaries@data[sel,"colour"], na.rm=TRUE), 3), sum(sel, na.rm=TRUE)), stringsAsFactors = FALSE)
    html <- kml_description(md, asText = TRUE, cwidth = 120, twidth = 240)
  }
  
  if(grid2poly == TRUE){
    pol <- grid2poly(obj@summaries["colour"])
  }

  kml_open(folder.name = folder.name, file.name = file.name)
  
  # add a description for the whole folder:
  kml.out <- get("kml.out", envir=plotKML.fileIO)
  description_txt <- sprintf('<description><![CDATA[%s]]></description>', html)
  parseXMLAndAdd(description_txt, parent=kml.out[["Document"]])  
  assign('kml.out', kml.out, envir=plotKML.fileIO)
  
  if(grid2poly == TRUE){  
    kml_layer(obj = pol, colour = colour, ...)
  }
  else {
    kml_layer(obj = obj@summaries, z.lim = c(0,1), colour = colour, raster_name = paste(folder.name, "_observed.png", sep=""), ...)
  }

  if(plot.svar==TRUE){
    kml_layer(obj = obj@summaries, colour = colour.sd, colour_scale = get("colour_scale_svar", envir = plotKML.opts), raster_name = paste(folder.name, "_observed.sd.png", sep=""), plot.legend = FALSE)  
  }
  
  # Realizations:
  for(i in 1:N.r){
    rel <- obj@realizations[[i]]
    kml_layer(obj = rel, TimeSpan.begin = i, TimeSpan.end = i+1)
  }

  # close the file:
  kml_close(file.name = file.name)
  if (kmz == TRUE){
      kml_compress(file.name = file.name)
  }
  # open KML file in the default browser:
  if(open.kml==TRUE){
    kml_View(file.name)
  } else {
    message(paste("Object written to:", file.name))
  }
    
})


# end of script;
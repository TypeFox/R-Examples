# Purpose        : Generic method to plot Species distribution models in Google Earth 
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Dev Status     : Alpha
# Note           : Only Google Earth 5.0 (and later) supports plain text content, as well as full HTML and JavaScript, within description balloons;

setMethod("plotKML", "SpatialMaxEntOutput", function(
  obj,
  folder.name = normalizeFilename(deparse(substitute(obj, env=parent.frame()))),
  file.name = paste(folder.name, ".kml", sep=""),
  html.file = obj@maxent@html,
  iframe.width = 800,
  iframe.height = 800,
  pngwidth = 280, 
  pngheight = 280,
  pngpointsize = 14,
  colour,
  shape = "http://plotkml.r-forge.r-project.org/icon17.png",
  kmz = get("kmz", envir = plotKML.opts),
  open.kml = TRUE,
  TimeSpan.begin = obj@TimeSpan.begin,
  TimeSpan.end = obj@TimeSpan.end,
  ...
){

  # Guess aesthetics if missing:
  if(missing(colour)){ 
    obj@sp.domain@data[,"colour"] <- obj@sp.domain@data[,1]
  } else {
    if(is.name(colour)|is.call(colour)){
      obj@sp.domain@data[,"colour"] <- eval(colour, obj@sp.domain@data)
    } else {
      obj@sp.domain@data[,"colour"] <- obj@sp.domain@data[,as.character(colour)]      
    }
  }

  # objects:
  spname <- obj@sciname
  pr <- as(obj@predicted[1], "SpatialPixelsDataFrame")
  names(pr) = "layer"

  # start writing the object:
  kml_open(folder.name = folder.name, file.name = file.name)
  kml.out <- get("kml.out", envir=plotKML.fileIO)
  
  # add a description for the whole folder:
  if(file.exists(html.file)){ html.file <- paste("file:///", html.file, sep="") }
  description_txt <- paste('<description><![CDATA[<iframe src="', html.file, '" width="', iframe.width, 'px" height="', iframe.height, 'px">]]></description>', sep="")
  parseXMLAndAdd(description_txt, parent=kml.out[["Document"]])  
  assign('kml.out', kml.out, envir=plotKML.fileIO)

  # occurrences:
  kml_layer.SpatialPoints(obj = obj@occurrences, shape = shape, TimeSpan.begin = TimeSpan.begin, TimeSpan.end = TimeSpan.end, labels = "", ...)
  # spatial domain (green colour):
  kml_layer(obj = obj@sp.domain, colour = colour, colour_scale = rgb(t(col2rgb(c("white", rep("dark green", 5))))/255), plot.legend = FALSE)
  # predicted values:
  kml_layer.SpatialPixels(obj = pr, colour = "layer", ...)

  # plot the contributions to the model:
  png(filename=paste(spname, "_var_contribution.png", sep=""), width=pngwidth, height=pngheight, bg="white", pointsize=pngpointsize)
  par(mar=c(4.5,4.5,.8,.8))
  .plot.maxent(obj@maxent, main = "")  
  dev.off()
  # add the plot:
  kml_screen(image.file = paste(spname, "_var_contribution.png", sep=""), position = "LL", sname = paste("Variable contribution for", spname))

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


# copied from the Dismo package (it does not export 'dismo::plot' method)
.plot.maxent <- function(x, sort=TRUE, main='Variable contribution', xlab='Percentage') {
		r <- x@results
		rnames <- rownames(r)
		i <- grep('.contribution', rnames)
		r <- r[i, ]
		names(r) <- gsub('.contribution', '', names(r))
		if (sort) {
			r <- sort(r)
		}
		dotchart(r, main=main, xlab=xlab)
}


# end of script;
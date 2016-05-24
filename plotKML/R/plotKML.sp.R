# Purpose        : Default methods to plot sp-type objects;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : ;
# Dev Status     : pre-Alpha
# Note           : these functions can be further customized;


setMethod("plotKML", "SpatialPointsDataFrame", function(obj, folder.name = normalizeFilename(deparse(substitute(obj, env=parent.frame()))), file.name = paste(folder.name, ".kml", sep=""), size, colour, points_names, shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png", metadata = NULL, kmz = get("kmz", envir = plotKML.opts), open.kml = TRUE, ...){

  # Guess aesthetics if missing:
  if(missing(size)){ 
    obj@data[,"size"] <- obj@data[,1]
  } else {
    if(is.name(size)|is.call(size)){
      obj@data[,"size"] <- eval(size, obj@data)
    } else {
      obj@data[,"size"] <- obj@data[,deparse(size)]      
    }
  }
  if(missing(colour)){ 
    obj@data[,"colour"] <- obj@data[,1] 
    message("Plotting the first variable on the list")  
  } else {
    if(is.name(colour)|is.call(colour)){
      obj@data[,"colour"] <- eval(colour, obj@data)
    } else {
      obj@data[,"colour"] <- obj@data[,as.character(colour)]      
    }
  }
  if(missing(points_names)){ 
    if(is.numeric(obj@data[,1])){ 
      points_names <- signif(obj@data[,1], 3) 
    } else {
      points_names <- paste(obj@data[,1])     
    }
  }
    
  # open for writing:
  kml_open(folder.name = folder.name, file.name = file.name)
 
  # write layer:
  if(is.numeric(obj@data[,"colour"])){
    kml_layer.SpatialPoints(obj, size = size, colour = colour, points_names = points_names, shape = shape, metadata = metadata, ...)
  } else {
    kml_layer.SpatialPoints(obj, colour = colour, points_names = points_names, shape = shape, metadata = metadata, ...)
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


setMethod("plotKML", "SpatialLinesDataFrame", function(obj, folder.name = normalizeFilename(deparse(substitute(obj, env=parent.frame()))), file.name = paste(folder.name, ".kml", sep=""), metadata = NULL, kmz = get("kmz", envir = plotKML.opts), open.kml = TRUE, ...){
   
  # open for writing:
  kml_open(folder.name = folder.name, file.name = file.name)
 
  # write layer:
  kml_layer.SpatialLines(obj, metadata = metadata, ...)

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


setMethod("plotKML", "SpatialPolygonsDataFrame", function(obj, folder.name = normalizeFilename(deparse(substitute(obj, env=parent.frame()))), file.name = paste(folder.name, ".kml", sep=""), colour, plot.labpt, labels, metadata = NULL, kmz = get("kmz", envir = plotKML.opts), open.kml = TRUE, ...){

  # Guess aesthetics if missing:
  if(missing(labels)){ 
    obj@data[,"labels"] <- obj@data[,1] 
  } else {
    if(is.name(labels)|is.call(labels)){
      obj@data[,"labels"] <- eval(labels, obj@data)
    } else {
      obj@data[,"labels"] <- obj@data[,deparse(labels)]      
    }
  }
  if(missing(colour)){ 
    obj@data[,"colour"] <- obj@data[,1]
    message("Plotting the first variable on the list")
  } else {
    if(is.name(colour)|is.call(colour)){
      obj@data[,"colour"] <- eval(colour, obj@data)
    } else {
      obj@data[,"colour"] <- obj@data[,as.character(colour)]      
    }
  }
  if(missing(plot.labpt)){ plot.labpt <- TRUE }
    
  # open for writing:
  kml_open(folder.name = folder.name, file.name = file.name)
 
  # write layer:
  kml_layer.SpatialPolygons(obj, colour = colour, plot.labpt = plot.labpt, labels = labels, metadata = metadata,  ...)

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

## Pixels Grids Raster
.plotKML.SpatialPixels <- function(obj, folder.name = normalizeFilename(deparse(substitute(obj, env=parent.frame()))), file.name = paste(folder.name, ".kml", sep=""), colour, raster_name, metadata = NULL, kmz = FALSE, open.kml = TRUE, ...){

  # the kml_layer.Raster works only with "Spatial" class:
  if(class(obj)=="RasterLayer"){
    obj <- as(obj, "SpatialGridDataFrame")
  }

  # Guess aesthetics if missing:
  if(missing(colour)){ 
    obj@data[,"colour"] <- obj@data[,1]
    message("Plotting the first variable on the list") 
  } else {
    if(is.name(colour)|is.call(colour)){
      obj@data[,"colour"] <- eval(colour, obj@data)
    } else {
      obj@data[,"colour"] <- obj@data[,as.character(colour)]      
    }
  }
  if(missing(raster_name)){ 
    raster_name <- paste(normalizeFilename(names(obj)[1]), ".png", sep="")
  }
   
  # open for writing:
  kml_open(folder.name = folder.name, file.name = file.name)
 
  # write layer:
  kml_layer(obj, colour = colour, raster_name = raster_name, metadata = metadata, ...)

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

}

setMethod("plotKML", "SpatialPixelsDataFrame", .plotKML.SpatialPixels)
setMethod("plotKML", "SpatialGridDataFrame", .plotKML.SpatialPixels)
setMethod("plotKML", "RasterLayer", .plotKML.SpatialPixels)


setMethod("plotKML", "SpatialPhotoOverlay", function(obj, folder.name = normalizeFilename(deparse(substitute(obj, env=parent.frame()))), file.name = paste(folder.name, ".kml", sep=""), dae.name, kmz = get("kmz", envir = plotKML.opts), open.kml = TRUE, ...){
  
  x <- strsplit(obj@filename, "/")[[1]]
  image.id <- x[length(x)]
  if(missing(dae.name)){ dae.name <- gsub(x=image.id, "\\.", "_") }
    
  # open for writing:
  kml_open(folder.name = folder.name, file.name = file.name, ...)
 
  # write layer:
  kml_layer.SpatialPhotoOverlay(obj, ...)

  # close the file:
  kml_close(file.name = file.name)
  if (kmz == TRUE){
      kml_compress(file.name = file.name, files=dae.name)
  }
  # open KML file in the default browser:
  if(open.kml==TRUE){
    kml_View(file.name)
  } else {
    message(paste("Object written to:", file.name))
  }

})


setMethod("plotKML", "SoilProfileCollection", function(obj, folder.name = normalizeFilename(deparse(substitute(obj, env=parent.frame()))), file.name = paste(folder.name, ".kml", sep=""), var.name, metadata = NULL, kmz = get("kmz", envir = plotKML.opts), open.kml = TRUE, ...){
  
  if(missing(var.name)){ var.name <- names(obj@horizons)[!(names(obj@horizons) %in% c(obj@idcol, obj@depthcols))][1] }
    
  # open for writing:
  kml_open(folder.name = folder.name, file.name = file.name)
 
  # write layer:
  kml_layer.SoilProfileCollection(obj, var.name = var.name, balloon = TRUE, metadata = metadata, ...)

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

## spacetime irregular vectors / spacetime full data frames...
.plotKML.ST <- function(obj, folder.name = normalizeFilename(deparse(substitute(obj, env=parent.frame()))), file.name = paste(folder.name, ".kml", sep=""), colour, shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png", points_names, kmz = get("kmz", envir = plotKML.opts), open.kml = TRUE, ...){

  # Guess aesthetics if missing:
  if(missing(colour)){ 
    obj@data[,"colour"] <- obj@data[,1]
    message("Plotting the first variable on the list") 
  } else {
    if(is.name(colour)|is.call(colour)){
      obj@data[,"colour"] <- eval(colour, obj@data)
    } else {
      obj@data[,"colour"] <- obj@data[,as.character(colour)]      
    }
  }
  if(missing(points_names)&class(obj@sp)=="SpatialPoints"){ 
    if(is.numeric(obj@data[,1])){ 
      points_names <- signif(obj@data[,1], 3) 
    } else {
      points_names <- paste(obj@data[,1])     
    }
  }
  # subset to target variable:
  obj <- obj[,,"colour"]
    
  # open for writing:
  kml_open(folder.name = folder.name, file.name = file.name)
 
  # write layer:
  if(class(obj@sp)=="SpatialPoints"){
    kml_layer(obj, shape = shape, colour = colour, points_names = points_names, ...)
  } else {
    kml_layer(obj, ...)  
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

}

setMethod("plotKML", "STIDF", .plotKML.ST)
setMethod("plotKML", "STFDF", function(obj, ...) .plotKML.ST(as(obj, "STIDF"), ...))
setMethod("plotKML", "STSDF", function(obj, ...) .plotKML.ST(as(obj, "STIDF"), ...))


## Trajectories:
setMethod("plotKML", "STTDF", function(obj, folder.name = normalizeFilename(deparse(substitute(obj, env=parent.frame()))), file.name = paste(folder.name, ".kml", sep=""), colour, start.icon = "http://maps.google.com/mapfiles/kml/pal2/icon18.png", kmz = get("kmz", envir = plotKML.opts), open.kml = TRUE, ...){
                            
  # Guess aesthetics if missing:
  if(missing(colour)){ 
    obj@data[,"colour"] <- obj@data[,1]
    message("Plotting the first variable on the list") 
  } else {
    if(is.name(colour)|is.call(colour)){
      obj@data[,"colour"] <- eval(colour, obj@data)
    } else {
      obj@data[,"colour"] <- obj@data[,as.character(colour)]      
    }
  }
    
  # open for writing:
  kml_open(folder.name = folder.name, file.name = file.name)
 
  # write layer:
  kml_layer.STTDF(obj, colour = colour, start.icon = start.icon, ...)

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


## List of objects of type "sp"
setMethod("plotKML", "list", function(obj, folder.name = normalizeFilename(deparse(substitute(obj, env=parent.frame()))), file.name = paste(folder.name, ".kml", sep=""), size = NULL, colour, points_names = "", shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png", plot.labpt = TRUE, labels = "", metadata = NULL, kmz = get("kmz", envir = plotKML.opts), open.kml = TRUE, ...){
   
  # check class of object:
  if(any(!(sapply(obj, class)=="SpatialPointsDataFrame"|sapply(obj, class)=="SpatialLinesDataFrame"|sapply(obj, class)=="SpatialPolygonsDataFrame"|sapply(obj, class)=="SpatialPixelsDataFrame"))){
    stop("List of objects of class SpatialPoints*, SpatialLines*, SpatialPolygons*, SpatialPixels* expected")
  }
    
  # open for writing:
  kml_open(folder.name = folder.name, file.name = file.name)

  ## target variable: 
  for(i in 1:length(obj)){    
  if(missing(colour)){ 
    obj[[i]]@data[,"colour"] <- obj[[i]]@data[,1] 
  } else {
    if(is.name(colour)|is.call(colour)){
      obj[[i]]@data[,"colour"] <- eval(colour, obj[[i]]@data)
    } else {
      obj[[i]]@data[,"colour"] <- obj[[i]]@data[,as.character(colour)]      
    }
    }
  }
    
  if(all(sapply(obj, class)=="SpatialPointsDataFrame")){
    for(i in 1:length(obj)){

      if(points_names == ""){ 
        if(is.numeric(obj[[i]]@data[,1])){ 
          points_names_i <- signif(obj[[i]]@data[,1], 3) 
        } else {
          points_names_i <- paste(obj[[i]]@data[,1])     
        }
      }

      if(is.numeric(obj[[i]]@data[,"colour"])){
        kml_layer.SpatialPoints(obj[[i]], colour = colour, points_names = points_names_i, shape = shape, metadata = metadata, ...)
      } else {
        kml_layer.SpatialPoints(obj[[i]], colour = colour, points_names = points_names, shape = shape, metadata = metadata, ...)
      }
    }
  }

  if(all(sapply(obj, class)=="SpatialLinesDataFrame")){    
    for(i in 1:length(obj)){
      kml_layer.SpatialLines(obj[[i]], metadata = metadata, ...)
    }
  }
    
  if(all(sapply(obj, class)=="SpatialPolygonsDataFrame")){
    for(i in 1:length(obj)){
      # Guess aesthetics if missing:
      if(labels == ""){ 
        labels_i <- obj[[i]]@data[,1] 
      } else {
        if(is.name(labels)|is.call(labels)){
          labels_i <- eval(labels, obj[[i]]@data)
        } else {
          labels_i <- obj[[i]]@data[,deparse(labels)]      
        }
      }
     
      kml_layer.SpatialPolygons(obj[[i]], colour = colour, plot.labpt = plot.labpt, labels = labels_i, metadata = metadata, ...)
    }
  }
    
  if(all(sapply(obj, class)=="SpatialPixelsDataFrame")){
    for(i in 1:length(obj)){
      
      ## subset to the bounding box if necessary:
      bbn <- round(diff(obj[[i]]@bbox[1,])/obj[[i]]@grid@cellsize[1])*round(diff(obj[[i]]@bbox[2,])/obj[[i]]@grid@cellsize[2])
      if(max(obj[[i]]@grid.index, na.rm=TRUE)>bbn){
        x <- as.data.frame(obj[[i]])
        suppressWarnings( gridded(x) <- ~x+y )
        proj4string(x) = obj[[i]]@proj4string
        obj[[i]] <- x
      }
      
      raster_name_i <- paste(names(obj[[i]])[1], "_", i, ".png", sep="")     
      kml_layer.SpatialPixels(obj[[i]], colour = colour, raster_name = raster_name_i, metadata = metadata, plot.legend=FALSE, ...)
    }
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
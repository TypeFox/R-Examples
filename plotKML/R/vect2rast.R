# Purpose        : Convert a vector map to a raster and (optional) write to a file;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : Dylan Beaudette (debeaudette@ucdavis.edu); Pierre Roudier (pierre.roudier@landcare.nz);
# Status         : tested
# Note           : The output pixel size is determined using simple cartographic principles (see [http://dx.doi.org/10.1016/j.cageo.2005.11.008]);


## POINTS
vect2rast.SpatialPoints <- function(obj, fname = names(obj)[1], cell.size, bbox, file.name, silent = FALSE, method = c("raster", "SAGA")[1], FIELD = 0, MULTIPLE = 1, LINE_TYPE = 0, GRID_TYPE = 2, ...){
  
    # add indicator value if missing data frame:
    if(!any(slotNames(obj) %in% "data")){
       obj <- SpatialPointsDataFrame(obj, data=data.frame(x=rep(1, length(obj))))
       names(obj) <- "mask"
    } else {
       obj <- obj[fname]
    }
    
    if(missing(bbox)) { bbox <- obj@bbox }
    # make an empty raster based on extent:
    if(missing(cell.size)){ 
    # print warning:
    if(length(obj)>10000){
    warning("Automated derivation of suitable cell size can be time consuming and can lead to artifacts.", immediate. = TRUE)
    }
    
      if(requireNamespace("spatstat", quietly = TRUE)){
        x <- as(obj, "ppp")
        nd <- spatstat::nndist(x$x, x$y)
        ndb <- boxplot(nd, plot=FALSE)
        cell.size <- signif(ndb$stats[3]/2, 2)
        if(cell.size==0){ stop("Estimated cell size is 0, consider removing duplicate points") }
      if(silent==FALSE){message(paste("Estimated nearest neighbour distance (point pattern):", cell.size*2))}
      }
    }

    if(method=="raster"){
      x <- GridTopology(cellcentre.offset=bbox[,1], cellsize=c(cell.size,cell.size), cells.dim=c(round(abs(diff(bbox[1,])/cell.size), 0), ncols=round(abs(diff(bbox[2,])/cell.size), 0)))
      r.sp <- SpatialGrid(x, proj4string = obj@proj4string)
      r <- raster(r.sp)
      # convert factors to integers:
      if(is.factor(obj@data[,1])){
       obj@data[,1] <- as.integer(obj@data[,1])
      }
      # rasterize - convert vector to a raster map:    
      obj <- obj[!is.na(obj@data[,1]),]
      in.r <- rasterize(obj, r, field = names(obj)[1], ...)
      res <- as(in.r, "SpatialGridDataFrame")
      names(res) = names(obj)[1]
      attr(res@bbox, "dimnames") = attr(obj@bbox, "dimnames")
      attr(res@grid@cellcentre.offset, "names") <- attr(obj@bbox, "dimnames")[[1]]
      attr(res@grid@cellsize, "names") <- attr(obj@bbox, "dimnames")[[1]]
      attr(res@grid@cells.dim, "names") <- attr(obj@bbox, "dimnames")[[1]]
    
      if(!missing(file.name)){
      writeRaster(in.r, filename=file.name, overwrite=TRUE)
      }
    }
    
    else{ 
    if(method=="SAGA"){   # SAGA GIS 2.0.8
   
      if(!rsaga.env()[["cmd"]]=="NULL"){
      
      tmf <- tempfile()
      tf <- set.file.extension(tmf, ".shp")
      if(requireNamespace("maptools", quietly = TRUE)){
        maptools::writePointsShape(obj, tf)
      } else {
        writeOGR(obj, tf, tmf, driver="ESRI Shapefile")
      }
      
      if(missing(file.name)){
        file.name <- set.file.extension(tempfile(), ".sgrd")
      }
      # rasterize map using SAGA GIS:
      rsaga.geoprocessor("grid_gridding", 0, param=list(INPUT=tf, FIELD=FIELD, MULTIPLE=MULTIPLE, LINE_TYPE=LINE_TYPE, GRID_TYPE=0, TARGET=0, USER_XMIN=bbox[1,1]+cell.size/2, USER_XMAX=bbox[1,2]-cell.size/2, USER_YMIN=bbox[2,1]+cell.size/2, USER_YMAX=bbox[2,2]-cell.size/2, USER_SIZE=cell.size, USER_GRID=file.name, GRID_TYPE=GRID_TYPE), show.output.on.console = silent)
      res <- readGDAL(set.file.extension(file.name, ".sdat"), silent = silent)
      names(res) = names(obj)[1]
      attr(res@bbox, "dimnames") = attr(obj@bbox, "dimnames")
      attr(res@grid@cellcentre.offset, "names") <- attr(obj@bbox, "dimnames")[[1]]
      attr(res@grid@cellsize, "names") <- attr(obj@bbox, "dimnames")[[1]]
      attr(res@grid@cells.dim, "names") <- attr(obj@bbox, "dimnames")[[1]]

      }
      else { stop("SAGA GIS path could not be located. See 'rsaga.env()' for more info.") }
    }  
    else{ stop("Method not available")
     }
    }
    
    return(res)
}

## LINES 
vect2rast.SpatialLines <- function(obj, fname = names(obj)[1], cell.size, bbox, file.name, silent = FALSE, method = c("raster", "SAGA")[1], FIELD = 0, MULTIPLE = 1, LINE_TYPE = 1, GRID_TYPE = 2, ...){

    # add indicator value if missing data frame:
    if(!any(slotNames(obj) %in% "data")){
       obj <- SpatialLinesDataFrame(obj, data=data.frame(x=rep(1, length(obj))))
       names(obj) <- "mask"
    } else {
       obj <- obj[fname]
    }
    
    if(missing(bbox)) { bbox <- obj@bbox }
    # make an empty raster based on extent:
    if(missing(cell.size)) { 
    # print warning:
    if(length(obj)>1000){
    warning("Automated derivation of suitable cell size can be time consuming and can lead to artifacts.", immediate. = TRUE)
    }
    
      if(requireNamespace("spatstat", quietly = TRUE)){
        x <- as(as(obj, "SpatialLines"), "psp")
        nd <- spatstat::nndist.psp(x)  # this can be time consuming!
        ndb <- boxplot(nd, plot=FALSE)
        cell.size <- signif(ndb$stats[3]/2, 2)
        if(cell.size==0){ stop("Estimated cell size is 0, consider removing duplicate points") }
      if(silent==FALSE){message(paste("Estimated nearest neighbour distance (line segments):", cell.size*2))}
      }
    }

    if(method=="raster"){
      x <- GridTopology(cellcentre.offset=bbox[,1], cellsize=c(cell.size,cell.size), cells.dim=c(round(abs(diff(bbox[1,])/cell.size), 0), ncols=round(abs(diff(bbox[2,])/cell.size), 0)))
      r.sp <- SpatialGrid(x, proj4string = obj@proj4string)
      r <- raster(r.sp)
      # convert factors to integers:
      if(is.factor(obj@data[,1])){
       obj@data[,1] <- as.integer(obj@data[,1])
      }
      # rasterize - convert vector to a raster map:    
      obj <- obj[!is.na(obj@data[,1]),]
      in.r <- rasterize(obj, r, field = names(obj)[1], ...)
      res <- as(in.r, "SpatialGridDataFrame")
      names(res) = names(obj)[1]
      attr(res@bbox, "dimnames") = attr(obj@bbox, "dimnames")
      attr(res@grid@cellcentre.offset, "names") <- attr(obj@bbox, "dimnames")[[1]]
      attr(res@grid@cellsize, "names") <- attr(obj@bbox, "dimnames")[[1]]
      attr(res@grid@cells.dim, "names") <- attr(obj@bbox, "dimnames")[[1]]  
    
      if(!missing(file.name)){
      writeRaster(in.r, filename=file.name, overwrite=TRUE)
      }
    }
    
    else{ 
    if(method=="SAGA"){   # SAGA GIS 2.0.8
   
      if(!rsaga.env()[["cmd"]]=="NULL"){
      
      tmf <- tempfile()
      tf <- set.file.extension(tmf, ".shp")
      if(requireNamespace("maptools", quietly = TRUE)){
        maptools::writeLinesShape(obj, tf)
      } else {
        writeOGR(obj, tf, tmf, driver="ESRI Shapefile")
      }

      if(missing(file.name)){
        file.name <- set.file.extension(tempfile(), ".sgrd")
      }
      # rasterize map using SAGA GIS:
      rsaga.geoprocessor("grid_gridding", 0, param=list(INPUT=tf, FIELD=FIELD, MULTIPLE=MULTIPLE, LINE_TYPE=LINE_TYPE, GRID_TYPE=0, TARGET=0, USER_XMIN=bbox[1,1]+cell.size/2, USER_XMAX=bbox[1,2]-cell.size/2, USER_YMIN=bbox[2,1]+cell.size/2, USER_YMAX=bbox[2,2]-cell.size/2, USER_SIZE=cell.size, USER_GRID=file.name, GRID_TYPE=GRID_TYPE), show.output.on.console = silent)
      res <- readGDAL(set.file.extension(file.name, ".sdat"), silent = silent)
      names(res) = names(obj)[1]
      attr(res@bbox, "dimnames") = attr(obj@bbox, "dimnames")
      attr(res@grid@cellcentre.offset, "names") <- attr(obj@bbox, "dimnames")[[1]]
      attr(res@grid@cellsize, "names") <- attr(obj@bbox, "dimnames")[[1]]
      attr(res@grid@cells.dim, "names") <- attr(obj@bbox, "dimnames")[[1]]
      
      }
      else { stop("SAGA GIS path could not be located. See 'rsaga.env()' for more info.") }
    }  
    else{ stop("Method not available")
     }
    }
    
    return(res)
}

## POLYGONS 
vect2rast.SpatialPolygons <- function(obj, fname = names(obj)[1], cell.size, bbox, file.name, silent = FALSE, method = c("raster", "SAGA")[1], FIELD = 0, MULTIPLE = 0, LINE_TYPE = 1, GRID_TYPE = 2, ...){
  
    # add indicator value if missing data frame:
    if(!any(slotNames(obj) %in% "data")){
       obj <- SpatialPolygonsDataFrame(obj, data=data.frame(x=rep(1, length(obj))))
       names(obj) <- "mask"
    } else {
       obj <- obj[fname]
    } 

    if(missing(bbox)) { bbox <- obj@bbox }
    # make an empty raster based on extent:
    if(missing(cell.size)) { 
    # print warning:
    if(length(obj)>1000){
    warning("Automated derivation of suitable cell size can be time consuming and can lead to artifacts.", immediate. = TRUE)
    }
    
      if(requireNamespace("spatstat", quietly = TRUE)){
        x <- sapply(obj@polygons, slot, "area")
        cell.size <- signif(sqrt(median(x))/2, 2)
        if(cell.size==0){ stop("Estimated cell size is 0, consider removing duplicate points") }
        if(silent==FALSE){message(paste("Estimated median polygon size:", cell.size*2))}
      }
    }

    if(method=="raster"){
      x <- GridTopology(cellcentre.offset=bbox[,1], cellsize=c(cell.size,cell.size), cells.dim=c(round(abs(diff(bbox[1,])/cell.size), 0), ncols=round(abs(diff(bbox[2,])/cell.size), 0)))
      r.sp <- SpatialGrid(x, proj4string = obj@proj4string)
      r <- raster(r.sp)
      # convert factors to integers:
      if(is.factor(obj@data[,1])){
       obj@data[,1] <- as.integer(obj@data[,1])
      }
      # rasterize - convert vector to a raster map:    
      obj <- obj[!is.na(obj@data[,1]),]
      in.r <- rasterize(obj, r, field = names(obj)[1], ...)
      res <- as(in.r, "SpatialGridDataFrame")
      names(res) = names(obj)[1]
      attr(res@bbox, "dimnames") = attr(obj@bbox, "dimnames")
      attr(res@grid@cellcentre.offset, "names") <- attr(obj@bbox, "dimnames")[[1]]
      attr(res@grid@cellsize, "names") <- attr(obj@bbox, "dimnames")[[1]]
      attr(res@grid@cells.dim, "names") <- attr(obj@bbox, "dimnames")[[1]]
    
      if(!missing(file.name)){
        writeRaster(in.r, filename=file.name, overwrite=TRUE)
      }
    }
    
    else{ 
    if(method=="SAGA"){   # SAGA GIS 2.0.8
   
      if(!rsaga.env()[["cmd"]]=="NULL"){
      
      tmf <- tempfile()
      tf <- set.file.extension(tmf, ".shp")
      if(requireNamespace("maptools", quietly = TRUE)){
        maptools::writePolyShape(obj, tf)
      } else {
        writeOGR(obj, tf, tmf, driver="ESRI Shapefile")
      }   
      
      if(missing(file.name)){
        file.name <- set.file.extension(tempfile(), ".sgrd")
      }
      # rasterize map using SAGA GIS:
      rsaga.geoprocessor("grid_gridding", 0, param=list(INPUT=tf, FIELD=FIELD, MULTIPLE=MULTIPLE, LINE_TYPE=LINE_TYPE, GRID_TYPE=0, TARGET=0, USER_XMIN=bbox[1,1]+cell.size/2, USER_XMAX=bbox[1,2]-cell.size/2, USER_YMIN=bbox[2,1]+cell.size/2, USER_YMAX=bbox[2,2]-cell.size/2, USER_SIZE=cell.size, USER_GRID=file.name, GRID_TYPE=GRID_TYPE), show.output.on.console = silent)
      res <- readGDAL(set.file.extension(file.name, ".sdat"), silent = silent)
      names(res) = names(obj)[1]
      attr(res@bbox, "dimnames") = attr(obj@bbox, "dimnames")
      attr(res@grid@cellcentre.offset, "names") <- attr(obj@bbox, "dimnames")[[1]]
      attr(res@grid@cellsize, "names") <- attr(obj@bbox, "dimnames")[[1]]
      attr(res@grid@cells.dim, "names") <- attr(obj@bbox, "dimnames")[[1]]
      
      }
      else { stop("SAGA GIS path could not be located. See 'rsaga.env()' for more info.") }
    }  
    else{ stop("Method not available")
     }
    }
    
    return(res)
}


setMethod("vect2rast", "SpatialPoints", vect2rast.SpatialPoints)
setMethod("vect2rast", "SpatialPolygons", vect2rast.SpatialPolygons)
setMethod("vect2rast", "SpatialLines", vect2rast.SpatialLines)

# end of script;

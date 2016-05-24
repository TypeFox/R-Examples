# Purpose        : Automatic reprojection of vector and raster features to geographic coordinates;
# Maintainer     : Pierre Roudier (pierre.roudier@landcare.nz);
# Contributions  : Tomislav Hengl (tom.hengl@wur.nl); Dylan Beaudette (debeaudette@ucdavis.edu); 
# Status         : tested
# Note           : in the case of gridded data, bounding box and cell size are estimated by the program (raster / GDAL);


reproject.SpatialPoints <- function(obj, CRS = get("ref_CRS", envir = plotKML.opts), ...) {
  message(paste("Reprojecting to", CRS, "..."))
  res <- spTransform(x = obj, CRSobj = CRS(CRS))
  return(res)
}


reproject.RasterLayer <- function(obj, CRS = get("ref_CRS", envir = plotKML.opts), program = "raster", tmp.file = TRUE, NAflag = get("NAflag", envir = plotKML.opts), show.output.on.console = FALSE, method, ...) {

  if(raster::is.factor(obj)){  
    method <- "ngb" 
  } else {  
    if(missing(method)){ method <- "bilinear" }
  }
  
  if(program=="raster"){
    message(paste("Reprojecting to", CRS, "..."))
    res <- raster::projectRaster(obj, crs = CRS, method = method, progress='text')
    names(res) <- names(obj)
  } else {
  
  if(program=="GDAL"){
    gdalwarp <- get("gdalwarp", envir = plotKML.opts)
  
    # look for GDAL path:  
    if(nchar(gdalwarp)==0){
      plotKML.env(silent = FALSE, show.env = FALSE)
      gdalwarp <- get("gdalwarp", envir = plotKML.opts)
    }
  
    if(tmp.file==TRUE){
       tf <- tempfile() 
    } else { 
       tf <- normalizeFilename(deparse(substitute(obj, env = parent.frame())))
    }
  
    if(!nchar(gdalwarp)==0){
      if(method == "ngb") { method <- "near" }
        writeRaster(obj, paste(tf, ".tif", sep=""), overwrite=TRUE, NAflag=NAflag)
        # resample to WGS84 system:
        message(paste("Using gdalwarp function:", gdalwarp))
        message(paste("Reprojecting to", CRS, "..."))
        system(paste(gdalwarp, " ", tf, ".tif", " -t_srs \"", CRS, "\" ", tf, "_ll.tif -dstnodata \"", NAflag, "\" ", " -r ", method, sep=""), show.output.on.console = show.output.on.console)
        res <- raster(paste(tf, "_ll.tif", sep=""), silent = TRUE)
        names(res) <- names(obj)
      } else {
        stop("Could not locate GDAL. See 'plotKML.env()' for more info.") }
  }
  }
  
  return(res)
}


reproject.SpatialGrid <- function(obj, CRS = get("ref_CRS", envir = plotKML.opts), tmp.file = TRUE, program = "raster", NAflag = get("NAflag", envir = plotKML.opts), show.output.on.console = FALSE, ...) {

  ## convert all character vectors to factors:
  for(j in 1:ncol(obj)){ 
    if(is.character(obj@data[,j])){ obj@data[,j] <- as.factor(obj@data[,j]) }
  }

  if(program=="raster"){

    # if multiple layers:
    if(ncol(obj) > 1) {
      r <- raster::stack(obj)
      res <- list(NULL)
      for(j in 1:ncol(obj)){
        if(is.factor(obj@data[,j])){
          r[[j]] <- raster::as.factor(r[[j]])
        }
        res[[j]] <- reproject(r[[j]], CRS = CRS) 
      }
      res <- as(raster::stack(res), "SpatialGridDataFrame")
      ## TH: time consuming but would be preferred:
      #  res <- as(res, "SpatialPixelsDataFrame")
      names(res) <- names(obj)
    }

    # single layer:
    else {
      r <- raster(obj)
      if(is.factor(obj@data[,1])){
        r <- raster::as.factor(r)
        message(paste("Reprojecting to", CRS, "..."))
        res <- as(raster::projectRaster(r, crs = CRS, method = "ngb"), "SpatialGridDataFrame")
      } else {
        res <- as(reproject(r, CRS = CRS), "SpatialGridDataFrame")
      }
      
      names(res) <- names(obj)
    }
    
    # try to fix factor-type objects:
    for(j in 1:ncol(obj)){
      if(is.factor(obj@data[,j])){
        # copy levels:
        res@data[,j] <- as.factor(res@data[,j])
        try( levels(res@data[,j]) <- levels(obj@data[,j]), silent = TRUE )
      }
    }
  }
  
  if(program=="GDAL"){
  gdalwarp <- get("gdalwarp", envir = plotKML.opts)
  # look for GDAL path if missing:  
  if(nchar(gdalwarp)==0){
    plotKML.env(silent = FALSE)
    gdalwarp <- get("gdalwarp", envir = plotKML.opts)
  }
  message(paste("Using gdalwarp function:", gdalwarp))
  
  if(!nchar(gdalwarp)==0){
  
    for(j in 1:ncol(obj)){
  
    if(tmp.file==TRUE){
        tf <- tempfile() 
        }
        else { 
          tf <- paste(normalizeFilename(deparse(substitute(obj, env = parent.frame()))), names(obj)[j], sep="_")
       }

        # write SPDF to a file:
        if(is.factor(obj@data[,j])){
          x <- obj[j]
          x@data[,1] <- as.integer(x@data[,1])
          if(max(x@data[,1], na.rm=TRUE)<254){
            writeGDAL(x, paste(tf, ".tif", sep=""), "GTiff", mvFlag = 255, type="Byte")
          } else {
            writeGDAL(x, paste(tf, ".tif", sep=""), "GTiff")
          }
        }        
        else {
          writeGDAL(obj[j], paste(tf, ".tif", sep=""), "GTiff")
        }
        
        message(paste("Reprojecting to", CRS, "..."))
        # resample to WGS84 system:
        if(is.factor(obj@data[,j])){
          system(paste(gdalwarp, ' ', tf, '.tif', ' -t_srs \"', CRS, '\" ', tf, '_ll.tif -dstnodata \"', NAflag, '\" -r near', sep=""), show.output.on.console = show.output.on.console)
        }
        else {
          system(paste(gdalwarp, ' ', tf, '.tif', ' -t_srs \"', CRS, '\" ', tf, '_ll.tif -dstnodata \"', NAflag, '\" -r bilinear', sep=""), show.output.on.console = show.output.on.console)
        }
        ## read images back to R:
        if(j==1){
          res <- readGDAL(paste(tf, "_ll.tif", sep=""), silent = TRUE)
          names(res) <- names(obj)[j]
        }
        else{
          res@data[names(obj)[j]] <- readGDAL(paste(tf, "_ll.tif", sep=""), silent = TRUE)$band1
        }
        
        ## reformat to the original factors:
          if(is.factor(obj@data[,j])){
            res@data[,j] <- as.factor(res@data[,j])
            try( levels(res@data[,j]) <- levels(obj@data[,j]) , silent = TRUE)
        }
        unlink(paste(tf, ".tif", sep=""))
        unlink(paste(tf, "_ll.tif", sep=""))        
      }
  } 
  
  else {
    stop("Could not locate GDAL. See 'plotKML.env()' for more info.") 
  }
  
  } 
  
  ## TH: time consuming but preferred:
  # res <- as(res, "SpatialPixelsDataFrame")
  return(res)
}


reproject.RasterStack <- function(obj, CRS = get("ref_CRS", envir = plotKML.opts)) {
  rs <- list(NULL)
  for(j in 1:nlayers(obj)){
    rs[[j]] <- reproject(obj[[j]], CRS = CRS)
  }
  return(stack(rs))
}


reproject.RasterBrick <- function(obj, CRS = get("ref_CRS", envir = plotKML.opts)) {
  r <- stack(obj)
  rs <- reproject.RasterStack(r, CRS = CRS)
  return(brick(rs))
}

# connect all methods and classes:
setMethod("reproject", "SpatialPoints", reproject.SpatialPoints)
setMethod("reproject", "SpatialPolygons", reproject.SpatialPoints)
setMethod("reproject", "SpatialLines", reproject.SpatialPoints)
setMethod("reproject", "RasterLayer", reproject.RasterLayer)
setMethod("reproject", "SpatialGridDataFrame", reproject.SpatialGrid)
setMethod("reproject", "SpatialPixelsDataFrame", reproject.SpatialGrid)
setMethod("reproject", "RasterStack", reproject.RasterStack)
setMethod("reproject", "RasterBrick", reproject.RasterBrick)

# end of script;

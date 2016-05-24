getSat <- function(path){
  band <- substr(list.files(path=path, pattern=paste0("^L[EC]\\d+\\w+\\d+_(B|band)2.(TIF|tif)$")), 0,3)
  if(length(band)==0){
    print(paste("ERROR: I expected something like landsat.band = LC82320832013319LGN00_BX.TIF in ", path))
    return()
  }
  if(band =="LC8"){return("L8")}
  if(band =="LE7"){return("L7")}
  if(band != "LE7" & band != "LC8"){
    print("Can't establish sat from band names")
    return()}
}


saveLoadClean <- function(imagestack, stack.names=NULL, file, ...){
  if(!getOption("waterWriteResults")){
    names(imagestack) <- stack.names
    return(imagestack)
  }
  tstamp <- NULL
  if(!getOption("waterOverwrite")){
    tstamp <- paste0("_", strftime(Sys.time(), format = "%Y%m%d%H%M%S"))}
  tmpdir <- getOption("waterOutputFolder")
  lastchar = substr(tmpdir, nchar(tmpdir), nchar(tmpdir))
  if (lastchar != "/" & lastchar != "\\") {
    tmpdir <- paste0(tmpdir, "/")
  }
  if(missing(file)){file <- paste0(deparse(substitute(imagestack)),".tif")}
  resultfile <- paste0(tmpdir, file, tstamp, ".tif")
  writeRaster(imagestack, filename = resultfile, ...)
  message(paste("Result saved as", paste0(tmpdir, file, tstamp, ".tif"), "    OK"))
  stack <- stack(resultfile)
  names(stack) <- stack.names
  removeTmpFiles(h=0)
  return(stack)
}


aoiCrop <- function(raster, aoi){
  if(!missing(aoi)){
    raster <- crop(raster,aoi)
    return(raster)
  }
  if(missing(aoi) & exists(x = "aoi", envir=.GlobalEnv)){
    aoi <- get(x = "aoi", envir=.GlobalEnv)
    raster <- crop(raster,aoi)
    return(raster)
  }
  return(raster)
}



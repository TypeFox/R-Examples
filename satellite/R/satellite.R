if ( !isGeneric("satellite") ) {
  setGeneric("satellite", function(x, ...)
    standardGeneric("satellite"))
}
#' Create a Satellite object
#'
#' @description
#' Method to create a Satellite object.
#' 
#' @param x A vector of filenames (see \code{raster::raster}) or a 
#' \code{raster::RasterStack}.
#' @param meta Optional metadata object (e.g. returned from 
#' \code{\link{compMetaLandsat}}). If x is a satellite dataset and recognised
#' as "Landsat", then the metadata is automatically extracted from the 
#' respective meta information file if both the satellite data and the metadata 
#' file follow the USGS Earth Explorer's naming convention.
#' @param log Optionally supply a log entry.
#' 
#' @return Satellite object
#' 
#' @export satellite
#' 
#' @details A satellite object consists of three data sections:
#' (i) a raster data section which holds the actual data values of the 
#' respective sensor bands, (ii) a metadata grid which holds meta information
#' for each sensor band (e.g. calibration coefficients, type of sensor band 
#' etc.) and (iii) a list of log information which records the processing 
#' history of the entire dataset.
#' 
#' @seealso \code{\link{compMetaLandsat}} to get more information about the
#' structure of the metadata component.
#' 
#' @name satellite
#' 
#' @examples
#' ## 'character' input (i.e. filenames)
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' 
#' satellite(files)
#' 
#' ## raster::RasterStack input
#' satellite(l8)
#' 
NULL


# Function using vector of filenames -------------------------------------------
#' 
#' @rdname satellite
#' 
setMethod("satellite", 
          signature(x = "character"), 
          function(x, meta, log){
            if(missing(meta)){
              if(lutInfoSGRPfromFilename(x) == "Landsat"){
                meta <- compMetaLandsat(x)
              } else {
                meta <- data.frame(LAYER = 
                                     tools::file_path_sans_ext(basename(x)),
                                   FILE = x)
              }
            }
            layers <- lapply(seq(nrow(meta)), function(y){
              act_layer <- raster(meta$FILE[y])
              names(act_layer) <- meta$BCDE[y]
              return(act_layer)
            })
            if(missing(log)){
              ps <- list(time = Sys.time(), info = "Initial import", 
                         layers = "all", output = "all")
              log <- list(ps0001 = ps)
            }
            #return(new("Satellite", layers = layers, meta = meta, log = log))
            tmp <- new("Satellite", layers = layers, meta = meta, log = log)
            return(addRasterMeta2Sat(tmp))
          })


# Function using readily existing raster layers --------------------------------
#'
#' @rdname satellite
#' 
setMethod("satellite", 
          signature(x = "RasterStack"), 
          function(x, meta, log){
            if(missing(meta)){
              meta <- data.frame(DATE = as.POSIXlt(Sys.Date(), tz = "UTC"),
                                 FILE = names(x))
            }
            # layers <- lapply(seq(nlayers(x)), function(y){
            #   x[[y]]
            # })
            layers <- raster::unstack(x)
            if(missing(log)){
              ps <- list(time = Sys.time(), info = "Initial import", 
                         layers = "all", output = "all")
              log <- list(ps0001 = ps)
            }
            #return(new("Satellite", layers = layers, meta = meta, log = log))
            tmp <- new("Satellite", layers = layers, meta = meta, log = log)
            return(addRasterMeta2Sat(tmp))
          })


# Function using list of raster layers -----------------------------------------
#' 
#' @rdname satellite
#' 
setMethod("satellite", 
          signature(x = "list"), 
          function(x, meta, log){
            if(missing(meta)){
              meta <- data.frame(DATE = as.POSIXlt(Sys.Date(), tz = "UTC"),
                                 FILE = names(x))
            }
            layers <- x
            if(missing(log)){
              ps <- list(time = Sys.time(), info = "Initial import", 
                         layers = "all", output = "all")
              log <- list(ps0001 = ps)
            }
            #return(new("Satellite", layers = layers, meta = meta, log = log))
            tmp <- new("Satellite", layers = layers, meta = meta, log = log)
            return(addRasterMeta2Sat(tmp))
          })

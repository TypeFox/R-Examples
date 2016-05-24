#' Get filename, bands and metadata file for Landsat 7 and 8 standard 1B/T format
#'
#' @description
#' The function compiles the sensor, band, filename and metadata filename information
#' for standard level 1B/T Landsat files.
#'
#' @param files Path and filename(s) of one or more Landsat band files or, 
#' alternatively, one or more Landsat metadata files. 
#'
#' @return \code{data.frame} containing filepaths, band numbers and metadata 
#' filepaths. 
#'
#' @export compFilePathLandsat
#'
#' @examples
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' compFilePathLandsat(files)  
#' 
compFilePathLandsat <- function(files){
  if((length(files) == 1 & grepl("MTL", files[1])) == FALSE){
    info <- lapply(files, function(x){
      layer <- tools::file_path_sans_ext(basename(x))
      pos <- gregexpr(pattern ='_B', layer)[[1]][1]
      band_ids <- substr(basename(x), pos + 2, 
                         nchar(layer))
      meta <- paste0(dirname(x), "/", substr(basename(x), 1, pos), "MTL.txt")
      sid <- substr(basename(x), 1, 3)
      
      sensor <- lutInfoSensorFromSID(sid)
      band_code <- lutInfoBCDEFromBID(sid = sid, bid = band_ids)
      data.frame(SID = sid, 
                 SENSOR = sensor,
                 BID = band_ids,
                 BCDE = band_code,
                 LAYER = layer,
                 FILE = x,
                 CALIB = "SC",
                 METAFILE = meta,
                 stringsAsFactors = FALSE,
                 row.names = NULL)
    })
    result <- (do.call("rbind", info))
    return(result)
  } else {
    info <- lapply(files, function(x){
      sid <- substr(basename(x), 1, 3)
      sensor <- lutInfoSensorFromSID(sid)
      band_code <- lutInfoBCDEFromBID(sid = sid)
      data.frame(BCDE = as.character(band_code),SID = sid, 
                 SENSOR = sensor,
                 BID = NA,
                 LAYER = NA,
                 FILE = NA,
                 CALIB = NA,
                 METAFILE = x,
                 stringsAsFactors = FALSE,
                 row.names = NULL)
    })
    result <- (do.call("rbind", info))
    rownames(result) <- NULL
    return(result)
  }
}
#' Get or access Satellite object information used by various functions
#'
#' @description
#' Get information from class Satellite.
#' 
#' @param data Data layer of a satellite object
#' @param type Type of the sensor band
#' @param lnbr Layer number
#' @param in_bcde BCDE of layer used as input dataset
#' @param info Log information added to metadata
#' @param meta_param Metadata parameters used to document new data layer
#' @param nbr Return specific data layer selected by number
#' @param out_bcde BCDE of layer used as output dataset
#' @param param Parameter of the metadata set (i.e. colname)
#' @param return_bcde Return bcde as attribute (TRUE/FALSE)
#' @param return_calib Return calibration information (TRUE/FALSE)
#' @param sat Satellite object (see \code{\link{satellite}}).
#' @param rst Input raster::Raster* object from which to extract metadata.
#' @param spectrum Spectral region, e.g. "solar" or "thermal".
#' @param width,flag Field width and format modifier for automated creation of 
#' BCDE information, defaults to '3' and '0', respectively. See 
#' \code{\link{formatC}} for further details.
#' @param prefix,postfix Prefix and postfix to be added to the created BCDE 
#' information.
#' @param calib Calibration information.
#' 
#' 
#' @return Objects of respective type (see \code{\link{satellite}}). 
#'
#' @details The functions are generally self-explaining in that sence that
#' \code{get*} returns the respective information and \code{set*} sets the
#' respective information from/in the Satellite object.
#'  
#' \code{addSatLog} adds a log entry to the Satellite object. 
#' 
#' @name satInfo
#' 
#' @examples
#' # List of input files
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#' 
#' # Raster stack l8
#' sat <- satellite(l8)
#' 
NULL


################################################################################
# Return (almost) complete sections of the Satellite object
################################################################################

#' @describeIn satInfo Return Satellite data layers
#' @export getSatDataLayers
#'
#'
getSatDataLayers <- function(sat, bcde = NULL){
  if (is.null(bcde)) {
    return(sat@layers)
  } else {
    ls_lyr <- sat@layers
    ch_lyr_nms <- sapply(ls_lyr, names)
    return(ls_lyr[ch_lyr_nms %in% bcde])
  }
}


# Return Satellite data layer i ------------------------------------------------
#' @export getSatDataLayer
#'
#' @describeIn satInfo Return Satellite data layer i
#'
getSatDataLayer <- function(sat, bcde){
  return(sat@layers[[getSatLNBR(sat, bcde)]])
}


# Return Satellite object metadata ---------------------------------------------
#' @export getSatMeta
#'
#' @describeIn satInfo Return Satellite object metadata
#'
getSatMeta <- function(sat, bcde){
  if(missing(bcde)){
    return(sat@meta)
  } else {
    return(sat@meta[sat@meta$BCDE %in% bcde, ])
  }
}


# Return template for Satellite object metadata which is based on existing band-
#' @export getSatMetaBCDETemplate
#'
#' @describeIn satInfo  Return template for Satellite object metadata which is based on existing band
#'
getSatMetaBCDETemplate <- function(sat, bcde){
  meta_template <- getSatMeta(sat, bcde)
  meta_template$DATE <- NULL
  meta_template$LAYER <- NULL
  meta_template$FILE <- NULL
  meta_template$METAFILE <- NULL
  return(meta_template)
}


# Return Satellite object log info ---------------------------------------------
#' @export getSatLog
#'
#' @describeIn satInfo Return Satellite object log info
#'
getSatLog <- function(sat){
  return(sat@log)
}

################################################################################
# Set entries to the Satellite object
################################################################################

# Set BCDE/data layer names of a Satellite object ------------------------------
#' @export setSatBCDE
#'
#' @describeIn satInfo Set BCDE/data layer names of a Satellite object
#'
setSatBCDE <- function(sat, bcde){
  ## if not supplied, BCDE is being created automatically
  if (missing(bcde))
    bcde <- createSatBCDE(sat)
  
  sat@meta$BCDE <- bcde
  for(i in seq(countSatDataLayers(sat))){
    names(sat@layers[[i]]) <- bcde[i]
  }
  return(sat)
}

# If not supplied, automatically create BCDE names of a Satellite object -------
#' @export createSatBCDE
#'
#' @describeIn satInfo If not supplied, automatically create BCDE names of a Satellite object
#'
createSatBCDE <- function(sat, width = 3, flag = 0, 
                          prefix = "B", postfix = "n") {
  int_nrow <- nrow(sat@meta)
  int_seq <- 1:int_nrow
  ch_seq <- formatC(int_seq, width = width, flag = flag)
  ch_bcde <- paste0(prefix, ch_seq, postfix)
  return(ch_bcde)
}


################################################################################
# Add entries to the Satellite object
################################################################################

# Add additional or overwrite metainformation parameter to Satellite object ----
#' @export addSatMetaParam
#'
#' @describeIn satInfo Add additional or overwrite metainformation parameter to Satellite object
#'
addSatMetaParam <- function(sat, meta_param){
  id <- colnames(meta_param)[1]
  name <- colnames(meta_param)[2]
  # Parameter already exists: overwrite, otherwise add
  if(length(which(name == colnames(sat@meta))) > 0){
    sat@meta[, which(name == colnames(sat@meta))] <- NULL
  } 
  sat@meta <- merge(sat@meta, meta_param, by = id, all.x = TRUE)
  sat@meta <- if (is.null(sat@meta$LNBR)) {
    sat@meta[order(sat@meta$BCDE), ]
    } else {
      sat@meta[order(sat@meta$LNBR), ]
    }
  return(sat)
}


# Add metainformation for an additional layer to Satellite object --------------
#' @export addSatMetaEntry
#'
#' @describeIn satInfo Add metainformation for an additional layer to Satellite object
#'
addSatMetaEntry <- function(sat, meta_param){
  if(missing(meta_param)){
    meta_param <- data.frame(DATE = as.POSIXlt(Sys.Date(), tz = "UTC"))
  }
  
  lnbr_new <- nrow(getSatMeta(sat)) + 1
  meta_param$LNBR <- lnbr_new
  
  if("DATE" %in% colnames(meta_param) == FALSE){
    meta_param$DATE <- as.POSIXlt(Sys.Date(), tz = "UTC")
  }
  
  if("LAYER" %in% colnames(meta_param) == FALSE){
    if(length(getSatDataLayers(sat)) == lnbr_new){
      meta_param$LAYER <- 
        getSatLayerfromData(sat, nbr = length(getSatDataLayers(sat)))
    } else {
      meta_param$LAYER <- paste0("Layer number ", lnbr_new)
    }
  }
  
  sat@meta <- plyr::rbind.fill(getSatMeta(sat), meta_param)
  return(sat)
}


# Add new log entry to Satellite object ----------------------------------------
#' @export addSatLog
#'
#' @describeIn satInfo Add new log entry to Satellite object
#'
addSatLog <- function(sat, info = NA_character_, in_bcde = NA_character_, 
                      out_bcde = NA_character_){
  new_length <- length(getSatLog(sat)) + 1
  ps <- sprintf("ps%04d", new_length)
  sat@log <- append(sat@log, list(list(time = Sys.time(), info = info, 
                                       in_bcde = in_bcde, out_bcde = out_bcde)))
  names(sat@log)[new_length] <- ps
  return(sat)
}


# Add new Satellite data layer --------------------------------------------------
#' @export addSatDataLayer
#'
#' @describeIn satInfo Add new Satellite data layer
#'
addSatDataLayer <- function(sat, bcde, data, meta_param, info, in_bcde){
  names(data) <- bcde
  sat@layers[[length(sat@layers) + 1]] <- data
  
  if(missing(meta_param)){
    meta_param = data.frame(BCDE = bcde)
  } else {
    meta_param$BCDE = bcde 
  }
  
  sat <- addSatMetaEntry(sat, meta_param = meta_param)
  sat <- addSatLog(sat, info = info, in_bcde = in_bcde, 
                   out_bcde = bcde)
  return(sat)
}


# Add raster meta data to Satellite object meta data ----
#' @export addRasterMeta2Sat
#'
#' @describeIn satInfo Add raster meta data to Satellite object meta data
#'
addRasterMeta2Sat <- function(sat){
  ## if BCDE is not available, it is automatically created 
  if (is.null(sat@meta$BCDE)) 
    sat <- setSatBCDE(sat)

    rst_meta <- data.frame(BCDE = sat@meta$BCDE, 
                         XRES = sapply(sat@layers, function(x) raster::xres(x)),
                         YRES = sapply(sat@layers, function(x) raster::yres(x)),
                         NROW = sapply(sat@layers, function(x) raster::nrow(x)),
                         NCOL = sapply(sat@layers, function(x) raster::ncol(x)),
                         NCELL = sapply(sat@layers, function(x) raster::ncell(x)),
                         XMIN = sapply(sat@layers, function(x) raster::xmin(x)),
                         XMAX = sapply(sat@layers, function(x) raster::xmax(x)),
                         YMIN = sapply(sat@layers, function(x) raster::ymin(x)),
                         YMAX = sapply(sat@layers, function(x) raster::ymax(x)),
                         PROJ = sapply(sat@layers, function(x) raster::projection(x)),
                         stringsAsFactors = FALSE)
  
  sat <- addSatMetaParam(sat, meta_param = rst_meta)
  
  
  return(sat)
}


# Create raster meta data ------------------------------------------------------
#' @export createRasterMetaData
#'
#' @describeIn satInfo Create raster meta data
#'
createRasterMetaData <- function(rst){
  rst_meta <- data.frame(XRES = raster::xres(rst),
                         YRES = raster::yres(rst),
                         NROW = raster::nrow(rst),
                         NCOL = raster::ncol(rst),
                         NCELL = raster::ncell(rst),
                         XMIN = raster::xmin(rst),
                         XMAX = raster::xmax(rst),
                         YMIN = raster::ymin(rst),
                         YMAX = raster::ymax(rst),
                         PROJ = raster::projection(rst),
                         stringsAsFactors = FALSE)
  
  return(rst_meta)
}


# Update raster meta data ------------------------------------------------------
#' @export updateRasterMetaData
#'
#' @describeIn satInfo Create raster meta data
#'
updateRasterMetaData <- function(sat, bcde) {
  
  rst_meta <- createRasterMetaData(getSatDataLayer(sat, bcde))
  
  for (i in colnames(rst_meta)) {
    sat@meta[sat@meta$BCDE == bcde, i] <- rst_meta[, i]
  }
   
  return(sat)
}



################################################################################
# Return individual entries of the Satellite object
################################################################################

# Return number of Satellite data layers ---------------------------------------
#' @export countSatDataLayers
#'
#' @describeIn satInfo Return number of Satellite data layers
#'
countSatDataLayers <- function(sat){
  return(length(sat@layers))
}


# Return parameter (general method implemented by the specific functions below)-
#' @param bcde Band code.
#' @export getSatParam
#'
#' @describeIn satInfo Return parameter (general method implemented by the specific functions below)
#' 
getSatParam <- function(sat, param, bcde, return_bcde = TRUE){
  if(length(which(param == colnames(getSatMeta(sat)))) > 0){
    if(param == "BCDE"){
      return(getSatMeta(sat)[, which(param == colnames(getSatMeta(sat)))])
    } else {
      if(missing(bcde)){
        param <- getSatMeta(sat)[, which(param == colnames(getSatMeta(sat)))]
        bcde <- as.character(getSatBCDE(sat))
      } else {
        param <- 
          getSatMeta(sat)[, 
                          which(param == colnames(getSatMeta(sat)))][match(
                            bcde, getSatMeta(sat)$BCDE)]
        bcde <- as.character(bcde)
      }
      if(return_bcde == TRUE){
        attr(param, "names") <- bcde
      }
      return(param)
    }
  } else {
    return(NA_character_)  
  }
}


# Return Band code -------------------------------------------------------------
#' 
#' @export getSatBCDE
#'
#' @describeIn satInfo Return Band code
#' 
getSatBCDE <- function(sat, lnbr){
  if(missing(lnbr)){
    getSatParam(sat, "BCDE", return_bcde = FALSE)  
  } else {
    getSatParam(sat, "BCDE", return_bcde = FALSE)[lnbr]
  }
}


# Return Band IDs --------------------------------------------------------------
#' 
#' @export getSatBID
#'
#' @describeIn satInfo Return Band IDs
#' 
getSatBID <- function(sat, bcde){
  getSatParam(sat, "BID", bcde, return_bcde = FALSE)
}


# Return sensor ID -------------------------------------------------------------
#' @export getSatSID
#'
#' @describeIn satInfo Return sensor ID
#' 
getSatSID <- function(sat){
  getSatParam(sat, "SID", return_bcde = FALSE)[1]
}


# Return sensor ----------------------------------------------------------------
#' @export getSatSensor
#'
#' @describeIn  satInfo Return sensor
#' 
getSatSensor <- function(sat){
  getSatParam(sat, "SENSOR", return_bcde = FALSE)[1]
}


# Return sensor group ----------------------------------------------------------
#' @export getSatSensorGroup
#'
#' @describeIn satInfo Return sensor group
#' 
getSatSensorGroup <- function(sat){
  getSatParam(sat, "SGRP", return_bcde = FALSE)[1]
}


# Return sensor information ----------------------------------------------------
#' @export getSatSensorInfo
#'
#' @describeIn satInfo Return sensor information
#' 
getSatSensorInfo <- function(sat){
  data.frame(SID = getSatSID(sat),
             SENSOR = getSatSensor(sat),
             SGRP = getSatSensorGroup(sat))
}


# Return spectrum --------------------------------------------------------------
#' @export getSatSpectrum
#'
#' @describeIn satInfo Return spectrum
#' 
getSatSpectrum <- function(sat, bcde){
  getSatParam(sat, "SPECTRUM", bcde)
}


# Return solar band codes ------------------------------------------------------
#' @export getSatBCDESolar
#'
#' @describeIn satInfo Return solar band codes
#' 
getSatBCDESolar <- function(sat){
  spectrum <- getSatSpectrum(sat)
  return(getSatBCDE(sat)[grep("solar", spectrum)])
}


# Return thermal band codes ------------------------------------------------------
#' @export getSatBCDEThermal
#'
#' @describeIn satInfo Return thermal band codes
#' 
getSatBCDEThermal <- function(sat){
  spectrum <- getSatSpectrum(sat)
  return(getSatBCDE(sat)[grep("thermal", spectrum)])
}


# Return sensor x resolution -----------------------------------------------------
#' @export getSatXRes
#'
#' @describeIn satInfo Return sensor x resolution
#' 
getSatXRes <- function(sat, bcde){
  getSatParam(sat, "XRES", bcde)
}


# Return sensor y resolution -----------------------------------------------------
#' @export getSatYRes
#'
#' @describeIn satInfo Return sensor y resolution
#' 
getSatYRes <- function(sat, bcde){
  getSatParam(sat, "YRES", bcde)
}


# Return mean sensor resolution (mean of x and y res) --------------------------
#' @export getSatRes
#'
#' @describeIn satInfo Return mean sensor resolution (mean of x and y res)
#' 
getSatRes <- function(sat, bcde){
  mean(getSatXRes(sat, bcde), getSatYRes(sat, bcde), na.rm = TRUE)
}


# Return sensor type -----------------------------------------------------------
#' @export getSatType
#'
#' @describeIn satInfo Return sensor type
#' 
getSatType <- function(sat, bcde){
  getSatParam(sat, "TYPE", bcde)
}


# Return CALIB -----------------------------------------------------------------
#' @export getSatCalib
#'
#' @describeIn satInfo Return calibration level
#' 
getSatCalib <- function(sat, bcde){
  getSatParam(sat, "CALIB", bcde)
}


# Return TYPE band codes matching id ------------------------------------------
#' @export getSatBCDEType
#'
#' @describeIn satInfo Return TYPE band codes
#' 
getSatBCDEType <- function(sat, bcde, type){
  type <- getSatType(sat, bcde)
  result <- getSatBCDE(sat)[grep(type, type)]
  if(length(result) == 0){
    result = NA_character_
  }
  return(result)
}


# Return BCDE matching TYPE ----------------------------------------------------
#' @export getSatBCDEFromType
#'
#' @describeIn satInfo Return BCDE matching TYPE
#' 
getSatBCDEFromType <- function(sat, type = "VIS"){
  as.character(stats::na.exclude(sat@meta$BCDE[sat@meta$TYPE == type]))
}


# Return BCDE matching TYPE ----------------------------------------------------
#' @export getSatBCDEFromSpectrum
#'
#' @describeIn satInfo Return BCDE matching TYPE
#' 
getSatBCDEFromSpectrum <- function(sat, spectrum = "solar"){
  as.character(stats::na.exclude(sat@meta$BCDE[sat@meta$SPECTRUM == spectrum]))
}



# Return SRES band codes matching type ------------------------------------------
#' @export getSatBCDESres
#'
#' @describeIn satInfo Return the mean of x and y resolution for band codes matching type
#' 
getSatBCDESres <- function(sat, bcde, type){
  sres <- getSatRes(sat, bcde)
  result <- getSatBCDE(sat)[grep(type, sres)]
  if(length(result) == 0){
    result = NA_character_
  }
  return(result)
}


# Return CALIB band codes matching type ------------------------------------------
#' @export getSatBCDECalib
#'
#' @describeIn satInfo Return calibration level for band codes matching type
#' 
getSatBCDECalib <- function(sat, bcde, calib){
  calib_val <- getSatCalib(sat, bcde)
  result <- getSatBCDE(sat)[grep(calib, calib_val)]
  if(length(result) == 0){
    result = NA_character_
  }
  return(result)
}


# Return CALIB band codes machting type and are solare bands ---------------------
#' @export getSatBCDESolarCalib
#'
#' @describeIn satInfo Return calibration level for band codes machting type and are solar bands
#' 
getSatBCDESolarCalib <- function(sat, bcde, calib){
  calib_val <- getSatBCDECalib(sat, bcde, calib)
  result <- getSatBCDESolar(sat)[getSatBCDESolar(sat) %in% calib_val]
  if(length(result) == 0){
    result = NA_character_
  }
  return(result)
}


# Return CALIB band codes machting type and are thermal bands --------------------
#' @export getSatBCDEThermalCalib
#'
#' @describeIn satInfo Return calibration level for band codes machting type and are thermal bands
#' 
getSatBCDEThermalCalib <- function(sat, bcde, calib){
  calib_val <- getSatBCDECalib(sat, bcde, calib)
  return(getSatBCDEThermal(sat)[getSatBCDEThermal(sat) %in% calib_val])
}


# Return band information ------------------------------------------------------
#' @export getSatBandInfo
#'
#' @describeIn satInfo Return band information
#' 
getSatBandInfo <- function(sat, bcde, return_calib = TRUE){
  if(return_calib){
    result <- data.frame(BID = getSatBID(sat, bcde),
                         #SRES = getSatRes(sat, bcde),
                         TYPE = getSatType(sat, bcde),
                         SPECTRUM = getSatSpectrum(sat, bcde),
                         CALIB = getSatCalib(sat, bcde))
  } else {
    result <- data.frame(BID = getSatBID(sat, bcde),
                         #SRES = getSatRes(sat, bcde),
                         TYPE = getSatType(sat, bcde),
                         SPECTRUM = getSatSpectrum(sat, bcde))
  }
  return(result)
}


# Return RAD_MAX ---------------------------------------------------------------
#' @export getSatRadMax
#'
#' @describeIn satInfo Return maximum radiance for bcde
#' 
getSatRadMax <- function(sat, bcde){
  getSatParam(sat, "RADMAX", bcde)
}


# Return RAD_MIN ---------------------------------------------------------------
#' @export getSatRadMin
#'
#' @describeIn satInfo Return minimum radiance for bcde
#' 
getSatRadMin <- function(sat, bcde){
  getSatParam(sat, "RADMIN", bcde)
}


# Return REF_MAX ---------------------------------------------------------------
#' @export getSatRefMax
#'
#' @describeIn satInfo Return maximum reflectance for bcde
#' 
getSatRefMax <- function(sat, bcde){
  getSatParam(sat, "REFMAX", bcde)
}


# Return REF_MIN ---------------------------------------------------------------
#' @export getSatRefMin
#'
#' @describeIn satInfo Return minimum reflectance for bcde
#' 
getSatRefMin <- function(sat, bcde){
  getSatParam(sat, "REFMIN", bcde)
}


# Return ESD -------------------------------------------------------------------
#' @export getSatESD
#'
#' @describeIn satInfo Return earth-sun distance
#' 
getSatESD <- function(sat){
  getSatParam(sat, "ESD")[1]
}


# Return ESun ------------------------------------------------------------------
#' @export getSatESUN
#'
#' @describeIn satInfo Return actual solar TOA irradiance
#' 
getSatESUN <- function(sat, bcde) {
  getSatParam(sat, "ESUN", bcde)
}


# Return SZEN ------------------------------------------------------------------
#' @export getSatSZEN
#'
#' @describeIn satInfo Return sun zenith angle
#' 
getSatSZEN <- function(sat, bcde){
  getSatParam(sat, "SZEN", bcde)
}


# Return SAZM ------------------------------------------------------------------
#' @export getSatSAZM
#'
#' @describeIn satInfo Return sun azimuth angle
#' 
getSatSAZM <- function(sat, bcde){
  getSatParam(sat, "SAZM", bcde)
}


# Return SELV ------------------------------------------------------------------
#' @export getSatSELV
#'
#' @describeIn satInfo Return Sun elevation
#' 
getSatSELV <- function(sat, bcde){
  getSatParam(sat, "SELV", bcde)
}

# Return Layer name from metadata ----------------------------------------------
#' @export getSatMetaLayer
#'
#' @describeIn satInfo Return Layer name from metadata
#' 
getSatMetaLayer <- function(sat, bcde){
  getSatParam(sat, "LAYER", bcde)
}


# Return Layer name from data layer --------------------------------------------
#' @export getSatLayerfromData
#'
#' @describeIn satInfo Return Layer name from data layer
#' 
getSatLayerfromData <- function(sat, bcde, nbr){
  if(missing(bcde)){
    layers <- sapply(getSatDataLayers(sat), function(x){
      names(x)
    })
    if(missing(nbr)){
      return(layers)
    } else {
      return(layers[nbr])
    }
  } else {
    names(getSatDataLayer(sat, bcde))
  }
}



# Return LNBR ------------------------------------------------------------------
#' @export getSatLNBR
#'
#' @describeIn satInfo Return Layer number
#' 
getSatLNBR <- function(sat, bcde){
  getSatParam(sat, "LNBR", bcde)
}


# Return LMin ------------------------------------------------------------------
#' @export getSatLMIN
#'
#' @describeIn satInfo Return minimum wavelength of the sensor band
#' 
getSatLMIN <- function(sat, bcde){
  getSatParam(sat, "LMIN", bcde)
}


# Return LMAX ------------------------------------------------------------------
#' @export getSatLMAX
#'
#' @describeIn satInfo Return maximum wavelength of the sensor band
#' 
getSatLMAX <- function(sat, bcde){
  getSatParam(sat, "LMAX", bcde)
}


# Return RADA ------------------------------------------------------------------
#' @export getSatRADA
#'
#' @describeIn satInfo Return addition coefficient for SC to radiance conversion
#' 
getSatRADA <- function(sat, bcde){
  getSatParam(sat, "RADA", bcde)
}


# Return RADM ------------------------------------------------------------------
#' @export getSatRADM
#'
#' @describeIn satInfo Return multiplicative coefficient for SC to radiance conversion
#' 
getSatRADM <- function(sat, bcde){
  getSatParam(sat, "RADM", bcde)
}


# Return REFA ------------------------------------------------------------------
#' @export getSatREFA
#'
#' @describeIn satInfo Return addition coefficient for SC to reflectance
#' 
getSatREFA <- function(sat, bcde){
  getSatParam(sat, "REFA", bcde)
}


# Return REFM ------------------------------------------------------------------
#' @export getSatREFM
#'
#' @describeIn satInfo Return multiplicative coefficient for SC to reflectance
#' 
getSatREFM <- function(sat, bcde){
  getSatParam(sat, "REFM", bcde)
}


# Return BTK1 ------------------------------------------------------------------
#' @export getSatBTK1
#'
#' @describeIn satInfo Return calibration coefficent to convert SC to brightness temperature
#' 
getSatBTK1 <- function(sat, bcde){
  getSatParam(sat, "BTK1", bcde)
}


# Return BTK2 ------------------------------------------------------------------
#' @export getSatBTK2
#'
#' @describeIn satInfo Return calibration coefficent to convert SC to brightness temperature
#' 
getSatBTK2 <- function(sat, bcde){
  getSatParam(sat, "BTK2", bcde)
}

# Return PRAD ------------------------------------------------------------------
#' @export getSatPRAD
#'
#' @rdname satInfo
#' 
getSatPRAD <- function(sat, bcde){
  getSatParam(sat, "PRAD", bcde)
}

# Return DATE ------------------------------------------------------------------
#' @export getSatDATE
#'
#' @describeIn satInfo Return DATE
#' 
getSatDATE <- function(sat, bcde){
  getSatParam(sat, "DATE", bcde)
}


# Return projection ------------------------------------------------------------------
#' @export getSatProjection
#'
#' @describeIn satInfo Return projection
#' 
getSatProjection <- function(sat, bcde){
  prj <- getSatParam(sat, "PROJ", bcde)
  
  if (is.na(prj))
    prj <- projection(sat@layers[[1]])
  
  return(prj)
}

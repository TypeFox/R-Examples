if (!isGeneric("convDN2RU") ) {
  setGeneric("convDN2RU", function(x, ...)
    standardGeneric("convDN2RU"))
}
#' Convert a band's scaled counts to radiance, reflectance and/or temperature
#'
#' @description
#' Convert a band's scaled counts to radiance, reflectance
#' and/or brightness temperature using a simple linear conversion without
#' any kind of atmospheric correction etc.
#' 
#' @param x An object of class Satellite, raster::RasterStack or 
#' raster::RasterLayer. 
#' @param convert Type of physical output; one of "rad", "ref", "bt" or "all". 
#' @param szen_correction Logical; if \code{TRUE}, sun zenith correction is 
#' being applied.
#' @param mult Multiplicative coefficient for value transformation (i.e. slope).
#' @param add Additive coefficient for value transformation (i.e. offset).
#' @param szen Cosine of solar zenith angle.
#' @param k1,k2 Temperature correction parameters.
#'   
#' @export convDN2RU
#' 
#' @name convDN2RU
#'
#' @details 
#' The conversion functions are taken from USGS' Landsat 8 manual
#' which is available online at 
#' \url{http://landsat.usgs.gov/Landsat8_Using_Product.php}.
#' 
#' @seealso \code{\link{calcAtmosCorr}} for conversions of scaled counts 
#' to physical units including a scene-based atmospheric correction.
#' 
#' @examples
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)  
#' sat <- convDN2RU(sat)
#' 
#' # If you use a raster layer, supply required meta information
#' bcde <- "B002n"
#' convDN2RU(x = getSatDataLayer(sat, bcde),
#'           mult = getSatRADM(sat, bcde),
#'           add = getSatRADA(sat, bcde))
#' 
NULL


# Function using satellite object ----------------------------------------------
#' 
#' @return If x is a Satellite object, a Satellite object with added converted 
#' layers; \cr 
#' if x is a \code{raster::Raster*} object, a \code{raster::Raster*} 
#' object with converted layer(s).
#' 
#' @rdname convDN2RU
#'
setMethod("convDN2RU", 
          signature(x = "Satellite"), 
          function(x, convert = "all", szen_correction = "TRUE"){
            .Deprecated("convSC2Rad", "convSC2Ref", "convRad2BT")
            if(convert == "all"){
              convert <- c("Rad", "Ref", "BT")
            }
            
            if("Rad" %in% convert){
              band_codes <- getSatBCDECalib(x, calib = "SC")
              for(bcde in band_codes){
                if(!is.na(getSatRADM(x, bcde))){
                  sensor_rad <- convDN2RU(x = getSatDataLayer(x, bcde),
                                              mult = getSatRADM(x, bcde),
                                              add = getSatRADA(x, bcde))
                  layer_bcde <- paste0(bcde, "_RAD")
                  
                  meta_param <- getSatMetaBCDETemplate(x, bcde)
                  meta_param$BCDE <- layer_bcde
                  meta_param$CALIB <- "RAD"
                  
                  info <- sys.calls()[[1]]
                  info <- paste0("Add layer from ", info[1], "(", 
                                 toString(info[2:length(info)]), ")")
                  
                  x <- addSatDataLayer(x, bcde = layer_bcde, data = sensor_rad,
                                       meta_param = meta_param,
                                       info = info, in_bcde = bcde)
                }
              }
            }
            
            if("Ref" %in% convert){
              band_codes <- getSatBCDESolarCalib(x, calib = "SC")
              for(bcde in band_codes){
                if(!is.na(getSatREFM(x, bcde))){
                  if(szen_correction == TRUE){
                    szen <- getSatSZEN(x, bcde)
                    sensor_ref <- convDN2RU(x = getSatDataLayer(x, bcde),
                                                mult = getSatREFM(x, bcde),
                                                add = getSatREFA(x, bcde),
                                                szen = szen)
                    calib = "REF"
                  } else {
                    sensor_ref <- convDN2RU(x = getSatDataLayer(x, bcde),
                                                mult = getSatREFM(x, bcde),
                                                add = getSatREFA(x, bcde))
                    calib = "REF_NoSZEN"
                  }
                  layer_bcde <- paste0(bcde, "_", calib)
                  
                  meta_param <- getSatMetaBCDETemplate(x, bcde)
                  meta_param$BCDE <- layer_bcde
                  meta_param$CALIB <- calib
                  
                  info <- sys.calls()[[1]]
                  info <- paste0("Add layer from ", info[1], "(", 
                                 toString(info[2:length(info)]), ")")
                  
                  x <- addSatDataLayer(x, bcde = layer_bcde, data = sensor_ref,
                                       meta_param = meta_param,
                                       info = info, in_bcde = bcde)
                }
              }
            }
            
            if("BT" %in% convert){
              band_codes <- getSatBCDEThermalCalib(x, calib = "SC")
              for(bcde in band_codes){
                if(!any(is.na(getSatRADM(x, bcde)), is.na(getSatBTK1(x, bcde)))){
                  sensor_ref <- convDN2RU(x = getSatDataLayer(x, bcde),
                                              mult = getSatRADM(x, bcde),
                                              add = getSatRADA(x, bcde),
                                              k1 = getSatBTK1(x, bcde),
                                              k2 = getSatBTK2(x, bcde))
                  layer_bcde <- paste0(bcde, "_BT")
                  
                  meta_param <- getSatMetaBCDETemplate(x, bcde)
                  meta_param$BCDE <- layer_bcde
                  meta_param$CALIB <- "BT"
                  
                  info <- sys.calls()[[1]]
                  info <- paste0("Add layer from ", info[1], "(", 
                                 toString(info[2:length(info)]), ")")
                  
                  x <- addSatDataLayer(x, bcde = layer_bcde, data = sensor_ref,
                                       meta_param = meta_param,
                                       info = info, in_bcde = bcde)
                }
              }
            }
            return(x)
          })


# Function using raster::RasterStack object ------------------------------------
#' 
#' @rdname convDN2RU
#'
setMethod("convDN2RU", 
          signature(x = "RasterStack"), 
          function(x, mult, add, szen, k1, k2){
            .Deprecated("convSC2Rad", "convSC2Ref", "convRad2BT")
            for(l in seq(nlayers(x))){
              x[[l]] <- convDN2RU(x[[l]], mult, add, szen, k1, k2)
            }
            return(x)
          })


# Function using raster::RasterLayer object ------------------------------------
#' 
#' @rdname convDN2RU
#'
setMethod("convDN2RU", 
          signature(x = "RasterLayer"), 
          function(x, mult, add, szen, k1, k2){
            .Deprecated("convSC2Rad", "convSC2Ref", "convRad2BT")
            result <- mult * x + add
            if(!missing(szen)){
              result <- result / cos(szen * pi / 180.0)
            }
            if(!missing(k1)){
              result <- k2 / log(k1 / result + 1)
            }
            return(result)
          })
if (!isGeneric("convRad2BT") ) {
  setGeneric("convRad2BT", function(x, ...)
    standardGeneric("convRad2BT"))
}
#' Convert a band's scaled counts to brightness temperature
#'
#' @description
#' Convert a band's radiance values to brightness temperature without
#' any kind of atmospheric correction etc.
#' 
#' @param x An object of class Satellite, raster::RasterStack or 
#' raster::RasterLayer providing radiance values.
#' @param k1,k2 Temperature correction parameters. 
#'   
#' @export convRad2BT
#' 
#' @name convRad2BT
#'
#' @details 
#' The conversion functions are taken from USGS' Landsat 8 manual
#' which is available online at 
#' \url{http://landsat.usgs.gov/Landsat8_Using_Product.php}.
#' 
#' @seealso \code{\link{calcAtmosCorr}} for converions of scaled counts 
#' to physical units including a scene-based atmospheric correction.
#' 
#' @examples
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)  
#' sat <- convRad2BT(sat)
#' 
NULL


# Function using satellite object ----------------------------------------------
#' 
#' @return If x is a Satellite object, a Satellite object with added converted 
#' layers; \cr
#' if x is a \code{raster::Raster*} object, a \code{raster::Raster*} object with
#' converted layer(s).
#' 
#' @rdname convRad2BT
#'
setMethod("convRad2BT", 
          signature(x = "Satellite"), 
          function(x){
            band_codes <- getSatBCDEThermalCalib(x, calib = "RAD")
            for(bcde in band_codes){
              if(!any(is.na(getSatRADM(x, bcde)), is.na(getSatBTK1(x, bcde)))){
                sensor_ref <- convRad2BT(x = getSatDataLayer(x, bcde),
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
            return(x)
          })


# Function using raster::RasterStack object ------------------------------------
#' 
#' @rdname convRad2BT
#'
setMethod("convRad2BT", 
          signature(x = "RasterStack"), 
          function(x, k1, k2){
            for(l in seq(nlayers(x))){
              x[[l]] <- convRad2BT(x[[l]], k1, k2)
            }
            return(x)
          })


# Function using raster::RasterLayer object ------------------------------------
#' 
#' @rdname convRad2BT
#'
setMethod("convRad2BT", 
          signature(x = "RasterLayer"), 
          function(x, k1, k2){
            x <- k2 / log(k1 / x + 1)
            return(x)
          })
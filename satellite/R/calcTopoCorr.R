if ( !isGeneric("calcTopoCorr") ) {
  setGeneric("calcTopoCorr", function(x, ...)
    standardGeneric("calcTopoCorr"))
}
#' Correct for topographic effects.
#' 
#' @param x Satellite object.
#' @param mask Logical. If \code{TRUE}, the cloudmask from the Satellite object 
#' (if available) will be considered in the regression model.
#' 
#' @details 
#' The method of Civco (1989) is applied on atmospherically corrected bands 
#' (if not already available in the Satellite object, 
#' \code{\link{calcAtmosCorr}} is performed with its default settings.):
#' First, an analytical hillshade image is created based on a DEM and sun 
#' elevation and sun zenith information from the metadata. A regression between
#' the hillshade (independent variable) and each channel is then calculated 
#' with consideration of a cloudmask (if available).
#' The regression coefficents are used to calibrate the hillshade raster 
#' (for each channel individually). 
#' Finally, the calibrated hillshade image is subtracted from the corresponding
#' channel and the mean value of the channel is added.
#' 
#' @return If x is a Satellite object, a Satellite object with added, 
#' topographic corrected layers; if x is a \code{raster::Raster*} object, a 
#' \code{raster::Raster*} object with converted layer(s).
#' 
#' @export calcTopoCorr
#' 
#' @name calcTopoCorr
#' 
#' @references CIVCO, D.L. (1989): Topographic normalization of Landsat Thematic
#' Mapper digitalimagery. \emph{Photogrammetric Engineering & Remote Sensing}, 
#' 55, 1303-1309.
#' 
#' @examples
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#' 
#' ## dem
#' files_dem <- list.files(path, pattern = "DEM", full.names = TRUE)
#' DEM <- raster(files_dem)
#' 
#' sat <- addSatDataLayer(sat, data = DEM, info = NULL, bcde = "DEM", in_bcde="DEM")
#' 
#' \dontrun{
#' sat <- calcTopoCorr(sat)
#' }
NULL


# Function using satellite object ----------------------------------------------
#' 
#' @rdname calcTopoCorr
#'
setMethod("calcTopoCorr", 
          signature(x = "Satellite"), 
          function(x, mask=TRUE){
            if (is.null(getSatDataLayer(x, "hillShade"))){
              if (is.null(getSatDataLayer(x, "DEM"))){stop(
                "please provide a DEM in the satellite object")}
              x <- demTools(x, method = "hillShade")
            }
            
            if(any(is.na(getSatBCDESolarCalib(x, calib = "AtmosCorr")))){
              x <- calcAtmosCorr(x)
            }
            
      #      atmoscbands <- getSatDataLayers(x)[grepl("_REF_AtmosCorr$", 
      #                                               getSatBCDE(x))]
            
            atmos_bands <- getSatBCDESolarCalib(x, calib = "AtmosCorr")
            
           
            for(i in 1:length(atmos_bands)){
              hillsh <- getSatDataLayer(x, "hillShade")
              #hillsh <- raster::resample (hillsh, atmoscbands[[i]])
              hillsh <- raster::resample (hillsh, getSatDataLayer(x, atmos_bands[i]))
              
              if (mask){
                cloudmask <- getSatDataLayer(x, "cloudmask")[[1]]
              }
              layer_bcde <- gsub("AtmosCorr","TopoCorr",getSatBCDE(x)[
                grepl("_REF_AtmosCorr$", getSatBCDE(x))][i])
              
              meta_param <- data.frame(getSatSensorInfo(x),
                                       getSatBandInfo(x, atmos_bands[i], 
                                                      return_calib = FALSE),
                                       CALIB = "REF_TopoCorr")
              
              tmp  <- calcTopoCorr(getSatDataLayer(x, atmos_bands[[i]]), hillsh, cloudmask)
              
              
              
              x <- addSatDataLayer(x, bcde = layer_bcde, data = tmp, 
                                   info="Add layer from calcTopoCorr(x)", 
                                   meta_param = meta_param,
                                   in_bcde=getSatBCDE(x)[
                                     grepl("_REF_AtmosCorr$", 
                                           getSatBCDE(x))][i])
            }
            return(x)
          })


# Function using raster::RasterStack object ------------------------------------
#' 
#' @rdname calcTopoCorr
#'
setMethod("calcTopoCorr", 
          signature(x = "RasterStack"), 
          function(x, hillsh, cloudmask=NULL){
            for(l in seq(nlayers(x))){
              x[[l]] <- calcTopoCorr(x[[l]], hillsh, cloudmask)
            }
            return(x)
          })



# Function using raster::RasterLayer object ------------------------------------
#' 
#' @param hillsh A \code{raster::RasterLayer} created with 
#' \code{raster::hillShade}. 
#' @param cloudmask A \code{raster::RasterLayer} in which clouds are masked with 
#' NA values. 
#' @rdname calcTopoCorr
#'
setMethod("calcTopoCorr", 
          signature(x = "RasterLayer"), 
          function(x, hillsh, cloudmask=NULL){
            xtmp <- x
            if (!is.null(cloudmask)){
              xtmp <- raster::mask(x, cloudmask)
            }
            model <- summary(stats::lm(raster::values(xtmp) ~ 
                                         raster::values(hillsh)))
            calib <- hillsh * model$coefficients[[2]] + model$coefficients[[1]]
            x <- x - calib + mean(raster::values(xtmp), na.rm=TRUE)
            return(x)
          })


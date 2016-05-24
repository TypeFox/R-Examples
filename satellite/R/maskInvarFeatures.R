if ( !isGeneric("maskInvarFeatures") ) {
  setGeneric("maskInvarFeatures", function(x, ...)
    standardGeneric("maskInvarFeatures"))
}
#' Identify pseudo-invariant features from a satellite scene
#'
#' @description
#' Identify pseudo-invariant features from a satellite scene based on a 
#' vis, near infravis and short-wave infravis band.
#'
#' @param x A Satellite object or a \code{raster::RasterLayer} providing the 
#' sensor's vis band.
#' @param nir A \code{raster::RasterLayer} containing the sensor's nir band.
#' @param swir A \code{raster::RasterLayer} containing the sensor's swir band.
#' @param quant A value v = [0...1] which is used to define the percentage
#' threshold values (thv) for invariant features (nir/vis ratio < thv, 
#' swir band values > 1-thv). 
#' 
#' @export maskInvarFeatures
#' 
#' @name maskInvarFeatures
#' 
#' @details Invariant features are identified as pixels which belong to the 
#' group of (i) the n lowest VIS/NIR ratios and of (ii) the highest n
#' SWIR values. The value of n is given by the parameter quant = [0...1].
#' 
#' @references This function is taken and only slightly modified from the PIF
#' function by Sarah C. Goslee (2011). Analyzing Remote Sensing Data in R: The 
#' landsat Package. Journal of Statistical Software,43(4), 1-25. URL 
#' \url{http://www.jstatsoft.org/v43/i04/}.
#' 
#' The underlying theory has been published by Schott RJ, Salvaggio C and 
#' Volchok WJ (1988) Radiometric scene normalization using pseudoinvariant 
#' features. Remote Sensing of Environment 26/1, 
#' doi:10.1016/0034-4257(88)90116-2, available online at
#' \url{http://www.sciencedirect.com/science/article/pii/0034425788901162}.
#'
#' @examples
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#' sat <- maskInvarFeatures(sat)
#' 
#' maskInvarFeatures(x = getSatDataLayer(sat, "B004n"), 
#'                   nir = getSatDataLayer(sat, "B005n"), 
#'                   swir = getSatDataLayer(sat, "B007n"))
#'
#' ## when dealing with a 'RasterStack'
#' rst <- stack(files[c(6, 7, 9)])
#' maskInvarFeatures(rst)
#' 
NULL


# Function using satellite object ----------------------------------------------
#' 
#' @return If x is a Satellite object, a Satellite object with added layer; \cr
#' if x is a \code{raster::RasterLayer} object, a a \code{raster::RasterLayer} 
#' object with added layers (1 indicates invariant pixels, 0 otherwise). 
#' 
#' @rdname maskInvarFeatures
#'
setMethod("maskInvarFeatures", 
          signature(x = "Satellite"), 
          function(x){
            bcde_vis <- "B004n"
            bcde_nir <- "B005n"
            bcde_swir <- "B007n"
            mask <- maskInvarFeatures(x = getSatDataLayer(x, bcde_vis), 
                                      nir = getSatDataLayer(x, bcde_nir), 
                                      swir = getSatDataLayer(x, bcde_swir))
            layer_bcde <- "M0000_InvarFeatures"
            
            meta_param <- getSatSensorInfo(x)
            meta_param$BCDE <- layer_bcde
            meta_param$CALIB <- "BINARY"
            
            info <- sys.calls()[[1]]
            info <- paste0("Add layer from ", info[1], "(", 
                           toString(info[2:length(info)]), ")")
            
            x <- addSatDataLayer(x, bcde = layer_bcde, data = mask,
                                 meta_param = meta_param,
                                 info = info, in_bcde = paste(bcde_vis, 
                                                              bcde_nir, 
                                                              bcde_swir,
                                                              sep = ", "))
            return(x)
          })



# Function using raster::RasterStack object ------------------------------------
#' 
#' @param id_vis Index of the visible band. 
#' @param id_nir Index of the near infravis band.
#' @param id_swir Index of the short-wave infravis band. 
#' 
#' @rdname maskInvarFeatures
#'
setMethod("maskInvarFeatures", 
          signature(x = "RasterStack"), 
          function(x, quant = 0.01, id_vis = 1L, id_nir = 2L, id_swir = 3L) {
            
            ## split 'RasterStack' into single layers
            vis <- x[[id_vis]]
            nir <- x[[id_nir]]
            swir <- x[[id_swir]]
            
            ratio_nir_vis <- nir/vis
            ratio_nir_vis_quant <- quantile(ratio_nir_vis, probs = quant, 
                                            na.rm = TRUE)
            swir_quant <- quantile(swir, probs = 1-quant, na.rm = TRUE)
            
            invar_feats <- ratio_nir_vis < ratio_nir_vis_quant & swir > swir_quant
            invar_feats[invar_feats == 0] <- NA
            return(invar_feats)
          })


# Function using raster::RasterLayer object ------------------------------------
#' 
#' @rdname maskInvarFeatures
#'
setMethod("maskInvarFeatures", 
          signature(x = "RasterLayer"), 
          function(x, nir, swir, quant=0.01) {
            ratio_nir_vis <- nir/x
            ratio_nir_vis_quant <- quantile(ratio_nir_vis, probs = quant, na.rm=TRUE)
            swir_quant <- quantile(swir, probs = 1-quant, na.rm=TRUE)
            
            invar_feats <- ratio_nir_vis < ratio_nir_vis_quant & swir > swir_quant
            invar_feats[invar_feats == 0] <- NA
            return(invar_feats)
          })

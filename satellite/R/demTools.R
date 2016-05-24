if ( !isGeneric("demTools") ) {
  setGeneric("demTools", function(x, ...)
    standardGeneric("demTools"))
}
#' Compute terrain characteristics from digital elevation models
#'
#' @description
#' Compute terrain characteristics from digital elevation models (DEM) using 
#' \code{raster::terrain} or \code{raster::hillShade}.
#' @param x A DEM provided as an object of class Satellite or RasterLayer.
#' @param method Currently "slope", "aspect" and "hillshade" are implemented.
#' @param bcde The name of the DEM layer in the satellite object. 
#' 
#' @seealso 
#' \code{raster::terrain}, \code{raster::hillShade}.
#' 
#' @export demTools
#' @name demTools
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
#' sat <- demTools(sat)
NULL




# Function using satellite object ----------------------------------------------
#' 
#' @return If x is a Satellite object, a Satellite object with added layer containing calculated 
#' terrain information; if x is a \code{raster::RasterLayer} object, a 
#' \code{raster::RasterLayer} object with calculated terrain information.
#' 
#' @rdname demTools
#'
setMethod("demTools", 
          signature(x  =  "Satellite"), 
          function(x,method = "hillShade",bcde = "DEM"){
            sunElev <- NULL
            sunAzim <- NULL
            if (method == "hillShade"){
              if (sum(stats::complete.cases(unique(getSatSELV(x)))) > 1 || 
                    sum(stats::complete.cases(unique(getSatSAZM(x))))>1){
                print("Warning: Satellite data have different Sun elevation or 
                      Sun azimuth values. Only the first element is used")}
              sunElev <- as.numeric(getSatSELV(x)[1])
              sunAzim <- as.numeric(getSatSAZM(x)[1])
            }
            result <- demTools(getSatDataLayer(x, bcde), sunElev = sunElev, 
                               sunAzim = sunAzim, 
                               method = method)
            x <- addSatDataLayer(x, bcde = method, data = result, info = paste0(
              "Add layer ", method), in_bcde = bcde)
          }
)




# Function using raster::RasterLayer object ------------------------------------
#' 
#' @param sunElev If \code{method = "hillShade"}, the elevation angle of the 
#' sun in degrees. See parameter \code{angle} in \code{\link{hillShade}}. 
#' @param sunAzim If \code{method = "hillShade"}, the sun azimuth angle in 
#' degree. See parameter \code{direction} in \code{\link{hillShade}}.
#' 
#' @rdname demTools
#'
setMethod("demTools", 
          signature(x  =  "RasterLayer"), 
          function(x, sunElev, sunAzim, method = "hillShade"){
            if(missing(sunElev) | missing(sunAzim))
              stop("Parameters 'sunElev' and/or 'sunAzim' missing.\n")
            if (method == "slope"){
              result <- raster::terrain(x, opt = "slope") 
            }
            if (method == "aspect"){
              result <- raster::terrain(x, opt = "aspect") 
            }
            if(method == "hillShade"){
              slope <- raster::terrain(x, opt  =  "slope")
              aspect <- raster::terrain(x, opt  =  "aspect")
              result <- raster::hillShade(slope = slope, aspect = aspect, 
                                  angle = sunElev, direction = sunAzim)
            }
            return(result)
          })
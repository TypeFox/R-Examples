if ( !isGeneric('projectSatellite') ) {
  setGeneric('projectSatellite', function(x, y, ...)
    standardGeneric('projectSatellite'))
}

#' Reproject a 'Satellite' object
#' 
#' @description
#' Reproject a satellite object. Either a \code{template} or \code{crs} 
#' must be supplied. If \code{crs} is not supplied, \code{\link{alignGeometry}}
#' is called.
#'
#' @param x Satellite or Raster* object to be resampled.
#' @param template Raster* or spatial data set from which geometry can be 
#' extracted.
#' @param type Type of bands (e.g. VIS, NIR) which should be considered. If not 
#' supplied, all types will be processed depending and bands to be processed can
#' be defined by band_codes.
#' @param band_codes Band ID(s) to be resampled. If not supplied and type is 
#' not given, too, all bands will be considered for resampling.
#' @param crs character or object of class 'CRS'. 
#' PROJ.4 description of the coordinate reference system. 
#' See \code{\link{projectRaster}} for details.
#' @param method Method for resampling; "bilinear" for bilinear interpolation 
#' (default) or "ngb" for nearest neighbor interpolation. See e.g. 
#' \code{\link{resample}}, \code{\link{projectRaster}}.
#'  
#' @export projectSatellite
#' 
#' @name projectSatellite
#' @aliases projectSatellite,Satellite-method
#' 
#' @examples
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#' 
#' projectSatellite(sat, crs = "+init=epsg:4326", band_codes = "B001n")

setMethod('projectSatellite', signature(x = 'Satellite',
                            y = 'ANY'), 
          function(x,
                   template, 
                   band_codes, 
                   type, 
                   crs,
                   method = c("bilinear", "ngb")) {
            
            if (missing(crs)) {
              alignGeometry(x, template, band_codes, type, method)
            } else {
              
              method <- method[1]
              
              if(!missing(type)){
                band_codes <- getSatBCDEFromType(x, type)
              }
              else if(missing(band_codes)) {
                band_codes <- getSatBCDE(x)
              } 
              
              for(bcde in band_codes) {
                pr <- projectRaster(from = getSatDataLayer(x, bcde),
                                    crs = crs, method = method)
               
                layer_bcde <- paste0(bcde, "_reprojected")
                meta_param <- getSatMetaBCDETemplate(x, bcde)
                meta_param$BCDE <- layer_bcde
#                 meta_param$XRES <- xres(pr)
#                 meta_param$YRES <- yres(pr)
                
                info <- sys.calls()[[1]]
                info <- paste0("Add layer from ", info[1], "(", 
                               toString(info[2:length(info)]), ")")
                x <- addSatDataLayer(x, bcde = layer_bcde, data = pr,
                                     meta_param = meta_param,
                                     info = info, in_bcde = bcde)
              }
              print(layer_bcde)
              x <- updateRasterMetaData(x, layer_bcde)
            }
            return(x)
          }
)

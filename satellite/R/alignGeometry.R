if ( !isGeneric("alignGeometry") ) {
  setGeneric("alignGeometry", function(x, ...)
    standardGeneric("alignGeometry"))
}
#' Align raster geometry between two data sets
#'
#' @description
#' Align raster data by bringing it in the same geometry and extent.
#' If the data set is not in the same projection as the template, the alignment
#' will be computed by reprojection. If the data has already the same
#' projection, the data set will be cropped and aggregated prior to resampling
#' in order to reduce computation time.
#'
#' @param x Satellite or Raster* object to be resampled.
#' @param template Raster* or spatial data set from which geometry can be 
#' extracted.
#' @param type Type of bands (e.g. VIS, NIR) which should be considered. If not 
#' supplied, all types will be processed depending and bands to be processed can
#' be defined by band_codes.
#' @param band_codes Band ID(s) to be resampled. If not supplied and type is 
#' not given, too, all bands will be considered for resampling.
#' @param method Method for resampling; "bilinear" for bilinear interpolation 
#' (default) or "ngb" for nearest neighbor interpolation. See e.g. 
#' \code{\link{resample}}, \code{\link{projectRaster}}.
#'  
#' @export alignGeometry
#' 
#' @name alignGeometry
#' @aliases alignGeometry,Satellite-method
#' 
#' @examples
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#' 
#' alignGeometry(sat, template = getSatDataLayer(sat, "B008n"), 
#'                band_codes = "B001n")


# Function using satellite object ----------------------------------------------
#' 
#' @return Satellite object with alligned geometries.
#' 
#' @rdname alignGeometry
#'
setMethod("alignGeometry", 
          signature(x = "Satellite"), 
          function(x, template, band_codes, type, method = c("bilinear", "ngb")){
            method <- method[1]
            if(!missing(type)){
              band_codes <- getSatBCDE(x)[which(getSatType(x) == type)]
            }
            else if(missing(band_codes)){
              band_codes <- getSatBCDE(x)
            } 
            for(bcde in band_codes){
              ag <- alignGeometry(x = getSatDataLayer(x, bcde),
                                  template = template, method = method)
              layer_bcde <- paste0(bcde, "_AG")
              meta_param <- getSatMetaBCDETemplate(x, bcde)
              meta_param$BCDE <- layer_bcde
              meta_param$XRES <- xres(template)
              meta_param$YRES <- yres(template)
              
              info <- sys.calls()[[1]]
              info <- paste0("Add layer from ", info[1], "(", 
                             toString(info[2:length(info)]), ")")
              x <- addSatDataLayer(x, bcde = layer_bcde, data = ag,
                                   meta_param = meta_param,
                                   info = info, in_bcde = bcde)
            }
            return(x)
          })


# Function using raster::RasterStack object ------------------------------------
#' 
#' @return raster::RasterStack object with alligned layers
#' 
#' @rdname alignGeometry
#'
setMethod("alignGeometry", 
          signature(x = "RasterStack"), 
          function(x, template, method = c("bilinear", "ngb")){
            method <- method[1]
            for(l in seq(nlayers(x))){
              x[[l]] <- alignGeometry(x[[l]], template, method)
            }
            return(x)
          })


# Function using raster::RasterLayer object ------------------------------------
#' 
#' @return raster::RasterLayer object with alligned layer
#' 
#' @rdname alignGeometry
#'
setMethod("alignGeometry", 
          signature(x = "RasterLayer"),
          function(x, template, method = c("bilinear", "ngb")){
            method <- method[1]
            if(raster::projection(x) == raster::projection(template)){
              x <- raster::crop(x, template, snap = "out")
              if(class(template) == "RasterLayer"){
                if(x@ncols / template@ncols >= 2){
                  factor <- floor(x@ncols/template@ncols)
                  x <- raster::aggregate(x, fact = factor, fun = mean, 
                                 expand = TRUE)
                }
                x <- raster::resample(x, template, method = method)
              }
            } else {
              x <- raster::projectRaster(x, template, method = method)
            }
            return(x)
          })

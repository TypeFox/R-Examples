if ( !isGeneric('stack') ) {
  setGeneric('stack', function(x, ...)
    standardGeneric('stack'))
}

#' convert selected layers of a Satellite object to a RasterStack
#' 
#' @description
#' Convert selected layers of a Satellite object to a RasterStack
#' 
#' @param x an object of class 'Satellite'
#' @param layer character vector (bcde codes) or integer vector (index) of 
#' the layers to be stacked
#' @param ... additional arguments passed on to \code{\link{stack}}
#' 
#' @examples 
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#' 
#' stck <- stack(sat, c("B001n", "B002n", "B003n"))
#' stck
#' 
#' @export
#' @name stack
#' @docType methods
#' @rdname stack
#' @aliases stack,Satellite-method
#' 

### RasterStack -----------------------------------------------------------
setMethod('stack', signature(x = 'Satellite'), 
          function(x,
                   layer = names(x),
                   ...) {
            
            if (is.numeric(layer)) layer <- names(x)[layer]
            
            lyrs <- getSatDataLayers(x, bcde = layer)
            xres <- getSatXRes(x, bcde = layer)
            if (length(unique(xres)) > 1) {
              skip <- which(xres == unique(xres)[2])
              warning("\nlayer ", names(skip), " has different resolution",
                      "\nnot stacking this layer")
            }
            ind <- duplicated(xres)
            ind[1] <- TRUE
            stck <- raster::stack(lyrs[ind], ...)
            
            return(stck)
          }
)


if ( !isGeneric("crop") ) {
  setGeneric("crop", function(x, ...)
    standardGeneric("crop"))
}

#' Crop satellite object
#'
#' @description
#' The function is a wrapper around the \code{\link{crop}} function to 
#' easily crop a Satellite object by an \code{\link{extent}} object.
#'
#' @param x Satellite object.
#' @param y \code{\link{extent}} object. 
#' @param subset Logical; if \code{TRUE} (default), all layers but the cropped 
#' ones are being dropped; if \code{FALSE}, cropped layers are appended to the 
#' Satellite object.
#'
#' @return A Satellite object consisting of cropped layers only. If 
#' \code{subset = FALSE}, a Satellite object with the cropped layers appended.
#' 
#' @export crop
#' 
#' @name crop
#' @aliases crop,Satellite-method
#'
#' @details Crop layers of a Satellite object to the size of a given 
#' \code{raster::extent} object.
#' 
#' @references Please refer to the respective functions for references.
#'  
#' @seealso This function is a wrapper for \code{raster::crop}.
#'
#' @examples
#' \dontrun{
#' ## sample data
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#'
#' ## geographic extent of georg-gassmann-stadium (utm 32-n)
#' ext_ggs <- raster::extent(484015, 484143, 5627835, 5628020)
#' 
#' ## crop satellite object by specified extent
#' sat_ggs <- crop(sat, ext_ggs)
#' 
#' plot(sat)
#' plot(sat_ggs)
#' }
setMethod("crop", 
          signature(x = "Satellite"), 
          function(x, y, subset = TRUE) {
            rad_bands <- getSatBCDE(x)
            for (bcde_rad in rad_bands) {
              ref <- crop(getSatDataLayer(x, bcde_rad), y)
              # keep all metadata except for file path since cropped 
              # layers are in memory and set calib column flag.
              meta_param <- getSatMeta(x, bcde_rad)
              meta_param$CALIB <- "cropped"
              meta_param$FILE <- NULL

              info <- sys.calls()[[1]]
              info <- paste0("Add layer from ", info[1], "(", 
                             toString(info[2:length(info)]), ")")
              x <- addSatDataLayer(x, bcde = bcde_rad, data = ref,
                                   meta_param = meta_param,
                                   info = info, in_bcde = bcde_rad)
            }
            
            if(subset == TRUE){
              x <- subset(x, cid = "cropped")
              #reset LNBR (dirty hack)
              x@meta$LNBR <- rep(1:nrow(x@meta))
              x@meta$CALIB <- "SC"
            }
            
            return(x)
          }
)

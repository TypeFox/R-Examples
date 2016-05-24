if ( !isGeneric("subset") ) {
  setGeneric("subset", function(x, ...)
    standardGeneric("subset"))
}
#' Subset Satellite object data layers
#'
#' @description
#' This function subsets a Satellite object and returns the extracted dataset
#' as Satellite object.
#'
#' @param x Satellite or raster::Raster* object providing the source band(s) to 
#' be adjusted.
#' @param sid Band numbers or bcde which should be extracted
#' @param cid Calibration information used for subsetting (only works if sid is
#' not supplied to the function)
#' @param i Layer index(es) for subsetting.
#'
#' @name subset
#' @export subset
#' @aliases subset,Satellite-method
#' 
#' @return A Satellite object
#' 
#' @examples
#' ## sample data
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#' 
#' sat[[2:5]]
#' subset(sat, cid = "SC")
NULL


# Function using satellite object ----------------------------------------------
#' 
#' @return A Satellite object
#' 
#' @rdname subset
#'
setMethod('subset', signature(x = 'Satellite'), 
          function(x, sid, cid) {
            if(!missing(sid)){
              if (!is.character(sid)) {
                sid <- getSatBCDE(x, sid)
              }
              data <- getSatDataLayers(x, sid)
              meta <- getSatMeta(x, sid)
              log <- getSatLog(x)
              x <- satellite(data, meta = meta, log = log)
            } else {
              if (is.character(cid)) {
                i <- subset(x@meta, x@meta$CALIB == cid)
                if (length(i)==0) {
                  stop('invalid layer names')
                } else if (length(i) < length(cid)) {
                  warning('invalid layer names omitted')
                }
                meta_cid <- i
                data_cid <- x@layers[as.integer(row.names(i))]
              } else {
                # cidting by row/list numbers makes only sense for multiples of 
                # channel numbers in case all channels are submitted to the 
                # satellite object in the first place. Therefore maybe checking 
                # for correct selection would need to be implemented for user 
                # friendlyness. Maybe defining sat object as list of obejcts, 
                # where each object is an instance of the sat object as it is 
                # now defined would make handling complete instances of sat 
                # objects on which some computation was applied easier to handle
                # (instead of appending them directly to the list and the meta 
                # data frame)?
                cid <- as.integer(cid)
                if (! all(cid %in% 1:length(x@layers))) {
                  stop('not a valid cid')
                }
                meta_cid <- x@meta[cid,]
                data_cid <- x@layers[cid]
              }
              x <- new("Satellite", layers = data_cid, meta = meta_cid)
            }
            return(x)
          })


# Function using satellite object ----------------------------------------------
#' 
#' @return A Satellite object
#' 
#' @rdname subset
#'
setMethod("[[", signature(x = "Satellite"), 
          function(x, i) {
            subset(x, i)
          }
)

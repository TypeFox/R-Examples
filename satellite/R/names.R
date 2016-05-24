if ( !isGeneric("names") ) {
  setGeneric("names", function(x, ...)
    standardGeneric("names"))
}
#' Get/set Satellite data layer names
#'
#' @description
#' Get/set Satellite data layer names, i.e. the BCDE id.
#' 
#' @param x A Satellite object.
#' @param value Band codes of the individual data layers.
#'
#' @export names
#' 
#' @name names
#' 
#' @examples
#' path <- system.file("extdata", package = "satellite")
#' files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
#' sat <- satellite(files)
#' names(sat)
#' new_names <- paste0(names(sat), "_test")
#' names(sat) <- new_names
NULL


# Function using satellite object ----------------------------------------------
#' 
#' @return Satellite data layer names as character vector.
#' 
#' @rdname names
#'
setMethod("names", 
          signature(x = "Satellite"), 
          function(x){
            as.character(getSatBCDE(x))
          }
)



# Function using satellite object ----------------------------------------------
#' 
#' @rdname names
#'
setMethod('names<-', 
          signature(x = "Satellite"), 
          function(x, value){
            nl <- countSatDataLayers(x)
            if (is.null(value)) {
              bcde <- rep('', nl)
            } else if (length(value) != nl) {
              stop('incorrect number of layer names')
            }
            x <- setSatBCDE(x, value)
            return(x)
          }
)
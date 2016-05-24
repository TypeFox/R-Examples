#'Function to check if a object is a gapdata class object
#'
#'This function checks if a object is a gapdata class object.
#'
#' @param x a object
#' @export is.gapdata
#' @aliases is.gapdata
#' @return a logical TRUE or FALSE
#' @keywords internal
#' 


#check if a gapdata
is.gapdata <- function(x){
  inherits(x, "gapdata")
}

#' Bertini Object Check
#'
#' Test whether an object is an bertini object.
#'
#' @param x object to be tested
#' @return Vector of logicals.
#' @export
#' @examples
#' 
#' # see ?bertini
#' 
#' 
#' 
is.bertini <- function(x){
  any(class(x) == "bertini")
}


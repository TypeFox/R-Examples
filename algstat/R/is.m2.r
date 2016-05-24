#' Macaulay2 Object Check
#'
#' Test whether an object is an m2 object.
#'
#' @param x object to be tested
#' @return Vector of logicals.
#' @export
#' @examples
#' \dontrun{
#' 
#' is.m2(m2("13^1000"))
#' 
#' }
#' 
is.m2 <- function(x){
  any(class(x) == "m2")
}


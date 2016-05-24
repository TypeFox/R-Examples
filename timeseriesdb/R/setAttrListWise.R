#' Set Attributes to Each Element of List According to a Given Vector
#'
#' An attribute is set to all elements of a list given a vector of possible instances of the
#' the attribute. Note that this function fails to excecute if the vector is not of the same length
#' list.
#'
#' @param li a list
#' @param attrib character name of the attribute 
#' @param vec vector containing all instances of the attribute
#' @export
setAttrListWise <- function(li,attrib,vec){
  stopifnot(length(li) == length(vec))
  for(i in seq_along(li)){
    attr(li[[i]],attrib) <- vec[i]
  }
  li
}

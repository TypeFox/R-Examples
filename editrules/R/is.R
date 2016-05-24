#' Check object class
#' @name is.editrules
#' @aliases is.editset is.editmatrix is.editarray
#' @param x object to be checked
#' @return \code{logical}
#' 
{}

#' @export
#' @rdname is.editrules
is.editset <- function(x) inherits(x,"editset")

#' @export
#' @rdname is.editrules
is.editmatrix <- function(x) return(inherits(x, "editmatrix"))

is.cateditmatrix <- function(x){
  return(inherits(x, "cateditmatrix"))
}

#' @export
#' @rdname is.editrules
is.editarray <- function(x) inherits(x,"editarray")






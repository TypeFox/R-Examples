#' @include editmatrix.R
#' @include editarray.R
#' @include editset.R
{}

#' convert to matrix
#'
#' @export
#' @method as.matrix editmatrix
#' @rdname editmatrix
#'
#' @return \code{as.matrix}: Augmented \code{matrix} of \code{editmatrix}. (See also \code{\link{getAb}}).
as.matrix.editmatrix <- function(x, ...){
   array(x, dim=dim(x), dimnames=dimnames(x))
}


#' convert to matrix
#' @export
#' @method as.matrix editarray
#' @rdname editarray
#' @return \code{as.matrix}: The boolean matrix part of the \code{editarray}.
as.matrix.editarray <- function(x,...){
    array(x,dim=dim(x),dimnames=dimnames(x))
}


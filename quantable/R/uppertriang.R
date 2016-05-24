#' get values of upper triangle from matrix
#' @export
#' @param mat matrix
#' @examples
#'
#' t = matrix(1:25,ncol=5)
#' uppertriang(t)
uppertriang <- function(mat){
  res<-mat[upper.tri(mat,diag=FALSE)]
  return( c(unlist(res)) )
}

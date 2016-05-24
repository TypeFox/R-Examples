#' Create an Indicator Matrix
#' 
#' Creates an indicator matrix from a grouping vector.
#' 
#' @param grp.vec Numeric vector giving the group membership.  
#' @param K Scalar indicating the number of groups. Defaults to the number of
#' unique elements in \code{grp.vec}.
#' @export indmat
indmat <- function(grp.vec, K = length(unique(grp.vec))){
  G <- diag(K)
  if(!is(grp.vec, 'numeric')) grp.vec <- as.numeric(grp.vec)
  G <- G[grp.vec, , drop = FALSE]
  colnames(G) <- paste("g",1:K,sep="")
  G
}

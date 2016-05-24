#' Calculate part of the RSS from one node to another.
#'
#' This is a helper function for RelationStrengthSimilarity that returns the component of RSS contributed by paths of one particular length r.
#' 
#' @param x numeric matrix, adjacency matrix where the [i,j] entry gives the strength of the link from node i to node j.
#' @param v1 numeric, index of the 'from' node.
#' @param v2 numeric, index of the 'to' node.
#' @param r numeric, length of paths examined from \code{v1} to \code{v2}.
#' @param prepped logical, whether or not the adjacency matrix \code{x} has had zeros entered on the diagonal and each row divided by the row sum.
#' @return numeric, the part of the Relation Strength Similarity score from \code{v1} to \code{v2} contributed by paths of length \code{r}.
#' @export
#' @seealso \code{\link{RelationStrengthSimilarity}}
#' @references
#' "Discovering Missing Links in Networks Using Similarity Measures", 
#' Hung-Hsuan Chen, Liang Gou, Xiaolong (Luke) Zhang, C. Lee Giles. 2012.
#' 
#' \url{https://github.com/shaptonstahl/}
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
#' @examples
#' M <- as.matrix(get.adjacency(graph.atlas(128)))
#' M
#' dils:::RssThisRadius(x=M, v1=5, v2=6, r=1)
#' dils:::RssThisRadius(x=M, v1=5, v2=6, r=2)
#' dils:::RssThisRadius(x=M, v1=5, v2=6, r=3)
#' dils:::RssThisRadius(x=M, v1=5, v2=6, r=4)
RssThisRadius <- function(x, v1, v2, r, prepped=FALSE) {
  if( FALSE == prepped ) {
    diag(x) <- 0
    x <- sweep(x, 1, rowSums(x), "/")
  }
  n <- nrow(x)
  
  if( v1 == v2 ) {
    out <- 0
  } else if( 1 == r ) {
    out <- x[v1, v2]
  } else if( 2 == r ) {
    out <- sum(x[v1,] * x[,v2])
  } else if( 3 == r) {
    y <- sapply(1:n, function(ell) {
      RssThisRadius(x, v1, ell, 2, prepped=TRUE) - x[v1, v2] * x[v2, ell]
    })
    out <- sum(x[,v2] * y) + x[v1, v2] * x[v2, v1] * x[v1, v2]
  } else if( 4 == r ) {
    y <- sapply(1:n, function(ell) {
      RssThisRadius(x, v1, ell, 3, prepped=TRUE) - 
        x[v2, ell] * sum(x[v1,] * x[,v2]) +
        x[v2, ell] * x[v1, ell] * x[ell, v2] +
        x[v1, v2] * x[v2, v1] * x[v1, ell] -
        x[v1, v2] * sum(x[v2,] * x[,ell])
    })
    out <- sum(x[,v2] * y) + 
      x[v1,v2] * x[v2,v1] * sum(x[v1,] * x[,v2]) +
      x[v1,v2] * x[v1,v2] * sum(x[v2,] * x[,v1])
  } else {
    stop("RssThisRadius not yet supported for this value of r")
  }
  return( out )
}

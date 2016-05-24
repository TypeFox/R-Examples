#' Calculate the RSS from one node to another.
#'
#' This is a helper function for RelationStrengthSimilarity that returns the RSS for a single directed dyad.
#' 
#' @param xadj numeric matrix, adjacency matrix where the [i,j] entry gives the strength of the link from node i to node j.
#' @param v1 numeric, index of the 'from' node.
#' @param v2 numeric, index of the 'to' node.
#' @param radius numeric, length of longest path examined from \code{v1} to \code{v2}.
#' @return numeric, the Relation Strength Similarity score from \code{v1} to \code{v2}.
#' @seealso \code{\link{RelationStrengthSimilarity}}
#' @references
#' "Discovering Missing Links in Networks Using Similarity Measures", 
#' Hung-Hsuan Chen, Liang Gou, Xiaolong (Luke) Zhang, C. Lee Giles. 2012.
#' 
#' \url{https://github.com/shaptonstahl/}
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
#' @details
#' This is an internal function.  There are no guardians and it assumes that
#' the adjacency matrix \code{xadj} has had zeros entered on the diagonal 
#' and then each row divided by the row mean.
#' @examples
#' M <- as.matrix(get.adjacency(graph.atlas(128)))
#' M
#' M <- sweep(M, 1, rowMeans(M), "/")
#' M
#' dils:::RssCell(xadj=M, v1=5, v2=6, radius=1)
#' dils:::RssCell(xadj=M, v1=5, v2=6, radius=2)
#' dils:::RssCell(xadj=M, v1=5, v2=6, radius=3)
#' dils:::RssCell(xadj=M, v1=5, v2=6, radius=4)
RssCell <- function(xadj, v1, v2, radius){
  out <- sum(sapply(1:radius, function(this.r) {
    RssThisRadius(xadj, v1, v2, this.r, prepped=TRUE)
  }))
  return( out )
}

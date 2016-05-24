#' Combine edgelists into a single data.frame
#'
#' Given two or more edgelists, create a single edgelist with multiple columns, two for the from and to nodes and one for the weights from each constituent network.
#' 
#' @param ... data.frames, edgelists to be merged.
#' @return data.frame, single multinetwork edgelist
#' @export
#' @seealso \code{\link{EdgelistFill}}
#' @references
#' \url{http://www.haptonstahl.org/R}
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
#' @examples
#' edgelist1 <- data.frame(expand.grid(letters[1:2], letters[1:2]), 
#'                         uniform=runif(4))
#' edgelist2 <- data.frame(v1=c("a", "a"), v2=c("a", "b"), manual=c(.3, .5))
#' MergeEdgelists(edgelist1, edgelist2)
MergeEdgelists <- function(...) {
  # Arguments to add: nodelist (provide manually and fill to this)
  #                   mergeNodeLists=FALSE (make nodelist the union of all the edgelists' nodelists
  #
  # retrieve ... arguments
  edgelists <- list(...)
  n.edgelists <- length(edgelists)
  
  # Guardians
  for(i in 1:n.edgelists) if(!is(edgelists[[i]], "data.frame")) stop("Each object to be merged must be a data.frame")
  if(0 == n.edgelists) stop("Two or more edgelists must be specified")
  if(1 == n.edgelists) return(edgelists[[1]])
  # Check for incompatibility
  nodelists <- list(sort(union(edgelists[[1]][,1], edgelists[[1]][,2])))
  for(i in 2:n.edgelists) {
    nodelists[[i]] <- sort(union(edgelists[[i]][,1], edgelists[[i]][,2]))
    if( !identical(nodelists[[1]], nodelists[[i]]) ) {
      stop("Edgelist ", i, " has a set of nodes that doesn't match those in edgelist 1.")
    }
  }
  
  # deal with default and missing values
  for(i in 1:n.edgelists) {
    edgelists[[i]] <- EdgelistFill(edgelists[[i]], nodelist=nodelists[[i]])
  }
  
  # perform the function
  out <- edgelists[[1]]
  for(i in 2:n.edgelists) {
    old.n.cols <- ncol(out)
    out <- cbind(out, edgelists[[1]][,-c(1:2)])
    names(out)[-c(1:old.n.cols)] <- names(edgelists[[i]])[-c(1:2)]
  }
  
  # prepare and return the output
  names(out) <- make.names(names(out), unique=TRUE)
  return(out)
}
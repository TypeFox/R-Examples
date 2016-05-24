#' Convert an adjacency matrix to filled edgelist.
#'
#' Given the adjacency matrix for a network return a data.frame listing all possible edges and the weights for each edge.
#' 
#' @param A matrix, see 'Details' for formatting assumptions.
#' @param nodelist character, optional list of node names.
#' @return data.frame, full list of all possible edges with weights for each in third column.
#' @export
#' @seealso \code{\link{EdgelistFromIgraph}}
#' @references
#' \url{https://github.com/shaptonstahl/}
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
#' @details This assumes that the row of the adjacency matrix indicates the node the edge is coming 'from', the column represent the node the edge is going 'to', and the value in the adjacency matrix is the weight given to the edge.
#' @examples
#' n <- 10
#' A <- matrix(rnorm(n*n), nrow=n)
#' A
#' EdgelistFromAdjacency(A)
#' 
#' n <- 100
#' A <- matrix(rnorm(n*n), nrow=n)
#' A
#' EdgelistFromAdjacency(A)
#' 
#' n <- 500
#' A <- matrix(rnorm(n*n), nrow=n)
#' A
#' \dontrun{EdgelistFromAdjacency(A)}
EdgelistFromAdjacency <- function(A, 
                                  nodelist=paste("node", 1:nrow(A), sep="")) {
  # Guardians
  stopifnot(nrow(A) == ncol(A),
            length(nodelist) == nrow(A),
            make.names(nodelist) == nodelist)
  
  n <- nrow(A)
  out <- expand.grid(nodelist, nodelist)
  out$weight <- 0
  names(out) <- c("fromnode", "tonode", "weight")
  for(j in 1:n) {
    out$weight[((j-1)*n+1):(j*n)] <- A[,j]
  }
  return(out)
}

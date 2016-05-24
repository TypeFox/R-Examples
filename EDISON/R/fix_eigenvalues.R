#' Modify network to ensure stationarity.
#' 
#' This function ensures that the eigenvalues of the network structure matrix
#' are smaller or equal to 1, thereby ensuring stationarity of the regression.
#' This is done by removing edges at random until the condition is satisfied.
#' 
#' 
#' @param network Original network structure, a matrix of size NumNodes by
#' NumNodes.
#' @param q Number of nodes.
#' @param gauss_weights If \code{TRUE}, use Gaussian regression weight, if
#' \code{FALSE} conserve original weights.
#' @return Returns a network with fewer eigenvalues than the original network,
#' but satisfying the stationarity condition.
#' @author Frank Dondelinger
#' @seealso \code{\link{generateNetwork}}
#' @export fix_eigenvalues
fix_eigenvalues <-
function(network, q, gauss_weights) {
  e_v = eigen(network);
  
  # Ignore the weights
  parents = abs(network) > 0;
  
  while(any(abs(e_v$values) > 1)) {
    
    edge_num = sum(parents); 
      
    weights_phase = matrix(rnorm(q*q, 0, 1), q, q);
    
    new_parents = parents;

    # Remove random edges until condition satisfied
    remove_changes = 1:edge_num %in% 
                       sample(1:edge_num, 1, replace=FALSE)
    non_edges = parents[parents > 0];
    non_edges[remove_changes] = 0;
    new_parents[parents>0] = non_edges;
      
    network = parents * weights_phase;
    
    e_v = eigen(network);
  }
  
  return(network);
}


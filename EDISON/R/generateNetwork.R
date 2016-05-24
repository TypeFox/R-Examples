#' Generate a random network.
#' 
#' This function generates a random network with changepoints for structure
#' changes, for simulating synthetic data.
#' 
#' 
#' @param lambda_2 Average number of parents for each node in the network
#' (parameter for a Poisson distribution).
#' @param q Number of nodes.
#' @param min_phase_length Minimum segment length.
#' @param k_bar Maximum number of changepoints. If \code{fixed=TRUE}, this is
#' equal to the number of changepoints.
#' @param l Length of the time series.
#' @param lambda_3 Average number of structure changes between two segments
#' (parameter for a Poisson distribution).
#' @param spacing \code{1} if segments are equally spaced, \code{0} if they are
#' spaced randomly (subject to the constraints of min_phase_length).
#' @param gauss_weights \code{1} if edge weights in the network are drawn from
#' N(0, 1), \code{0} if they are fixed to be 1.
#' @param same \code{1} if all segments have the same network structure (no
#' changes), \code{0} otherwise.
#' @param change_method \code{'sequential'} if the changes happen sequentially
#' (i.e. changes at segment i are applied to segment i-1),
#' \code{'hierarchical'} if the changes happen with respect to a hypernetwork
#' (i.e. changes at segment i are applied to segment 0).
#' @param fixed \code{T} if the changepoint locations are fixed, \code{F} if
#' they should be sampled.
#' @param cps Changepoint locations (if they are fixed).
#' @return A list with the following elements: \item{network}{The network, a
#' list of length NumSegs, where each element is a NumNodes by NumNodes
#' matrix.} \item{epsilon}{The vector of changepoint locations.} \item{k}{The
#' number of changepoint.} \item{changes}{The number of changes among
#' segments.}
#' @author Frank Dondelinger
#' @seealso \code{\link{simulateNetwork}}
#' @examples
#' 
#' # Generate random network with default parameters
#' network = generateNetwork()
#' 
#' # Simulate data using generated network
#' dataset = simulateNetwork(net=network)
#' 
#' # Generate random network with 4 changepoints and 15 nodes, 
#' # with changepoints distributed over a timeseries of length 50
#' network = generateNetwork(l=50, q=15, fixed=TRUE, k_bar=4)
#' 
#' # Simulate data of length 50 using generated network
#' dataset = simulateNetwork(net=network)
#' 
#' @export generateNetwork
generateNetwork <-
function(lambda_2=0.45, q=10, min_phase_length=1, k_bar=5, l=10, 
                            lambda_3=2, spacing=1, gauss_weights=TRUE, same=FALSE, 
                            change_method='sequential',
                            fixed=FALSE, cps=NULL) { 
  # Draw hyperparameters
  lambda_1 = 6  

  k = k_bar + 1;
  epsilon = c();

  if(!is.null(cps)) {
    epsilon = cps;
    k = length(epsilon);
  } else if(!fixed) {
    # Choose number and location of change points
    while(k > k_bar || any(c(epsilon, l) - c(0, epsilon) < min_phase_length)) {
      k = rpois(1, lambda_1);
    
      if(spacing == 1) {
    	  epsilon = round(seq(round(l/k), l+1, length.out=(k+1)))
    	  epsilon = epsilon[-length(epsilon)]
      } else {
        epsilon = sort(sample(1:l, k, replace=FALSE))
      }
    }
  } else {
    k = k - 1
    if(spacing == 1) {
      epsilon = round(seq(round(l/(k+1)), l+1, length.out=(k+1)))
      epsilon = epsilon[-length(epsilon)]
    } else {
      epsilon = sort(sample(1:l, k, replace=FALSE))
    
      while(any(c(epsilon, l) - c(0, epsilon) < min_phase_length)) {
        epsilon = sort(sample(1:l, k, replace=FALSE))
      }
    }
  }
  
  # Draw initial network from prior
  parents = matrix(FALSE, q, q);

  for(i in 1:q) {

    num_parents = q+1;

    while(num_parents > q) {
      num_parents = rpois(1, lambda_2); 
    }
  
    parents_i = 1:q %in% sample(1:q, num_parents, replace=FALSE) 
	  parents[i, parents_i] = TRUE; 
    
    parents[i, i] = TRUE;
  }

  parents = t(parents);

  # Make sure original network generates stable data
  weights_orig = matrix(rnorm(q*q, 0, 1), q, q);
  
  network = fix_eigenvalues(parents*weights_orig, q, gauss_weights)
  
  parents = abs(network) > 0;
  
  # Original network (think of this as segment 0)
  parents_orig = parents;
  
  parents_prev = parents;

  network = list(); changes = list()

  # For each phase, draw number of changes and changed edges with 
  # respect to previous phase, then record new network  
  for(i in 1:(k+1)) {  
    
    # For a hierarchical model, apply changes with respect to the original 
    # network, rather than the previous segment.
    if(change_method == 'hierarchical') {
      parents = parents_orig
      parents_prev = parents_orig
    }
    
    for(j in 1:q) {
      
      change_num = rpois(1, lambda_3);
     
      # Make sure changes are applied equally to edges and non-edges 
      if(change_num > 0.5*q)
        change_num = round(0.5*q)
    
      add_edges = sum(runif(change_num) > 0.6);
      remove_edges = change_num - add_edges;      
      parents_j = parents[,j];
      edge_num = sum(parents_j);
      non_edge_num = q - edge_num;  
    
      
      if(remove_edges > edge_num)
        remove_edges = edge_num
        
      if(edge_num > 0.5*q)
        add_edges = 0;
        
      new_parents_j = parents_j;
    
      if(add_edges > 0) {
        add_changes = 1:non_edge_num %in% sample(1:non_edge_num, add_edges, replace=FALSE)
        edges = parents_j[parents_j == 0];
        edges[add_changes] = 1;
        new_parents_j[parents_j == 0] = edges;
      }
    
      if(remove_edges > 0) {
        remove_changes = 1:edge_num %in% sample(1:edge_num, remove_edges, replace=FALSE)
        non_edges = parents_j[parents_j > 0];
        non_edges[remove_changes] = 0;
        new_parents_j[parents_j>0] = non_edges;
      }
      
      if(!same) {
        parents[,j] = new_parents_j;
      }
    }
    
    parent_num = colSums(parents)
    parent_num[parent_num == 0] = 1;
    
    weights_phase = matrix(rnorm(q*q, 0, 1), q, q);
   
    # Apply weights if needed 
    if(gauss_weights) {
      network[[i]] = parents * weights_phase;
    } else {
      network[[i]] = parents 
    }
    
    parents = abs(network[[i]]) > 0;

    changes[[i]] = sum(abs(parents_prev - parents))
    parents_prev = parents
    
  }
  
  epsilon = c(epsilon, l);
 
  return(list(network=network, epsilon=epsilon, k=k, changes=changes, l=l))
}


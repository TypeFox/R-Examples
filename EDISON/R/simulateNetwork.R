#' Generate network and simulate data.
#' 
#' This function generates a random network with structure changepoints (or
#' takes one as input) and simulated data from it using a regression model.
#' 
#' 
#' @param l Length of the time series.
#' @param min_phase_length Minimum segment length.
#' @param k_bar Maximum number of changepoints.
#' @param q Number of nodes.
#' @param lambda_2 Average number of parents for each node in the network
#' (parameter for a Poisson distribution).
#' @param noise Standard deviation of the Gaussian observation noise. Can be
#' constant, or segment specific (in which case the number of changepoints
#' needs to be fixed and the noise needs to be a vector of the same length).
#' @param net Input network, can be \code{NULL} if a new network should be
#' generated.
#' @param lambda_3 Average number of structure changes between two segments
#' (parameter for a Poisson distribution).
#' @param spacing \code{1} if segments are equally spaced, \code{0} if they are
#' spaced randomly (subject to the constraints of min_phase_length).
#' @param gauss_weights \code{1} if edge weights in the network are drawn from
#' N(0, 1), \code{0} if they are fixed to be 1.
#' @param same \code{1} if the networks should all be the same (no changes),
#' \code{0} otherwise.
#' @param changes \code{'sequential'} if the changes happen sequentially (i.e.
#' changes at segment i are applied to segment i-1), \code{'hierarchical'} if
#' the changes happen with respect to a hypernetwork (i.e. changes at segment i
#' are applied to segment 0).
#' @param fixed \code{T} if the changepoint locations are fixed, \code{F} if
#' they should be sampled.
#' @param cps Changepoint locations (if they are fixed).
#' @param saveFile If not \code{NULL}, indicates the filename for saving the
#' output in R data format.
#' @return A list with elements: \item{sim_data}{A matrix of length NumNodes by
#' NumTimepoints containing the simulated data from the regression model.}
#' \item{epsilon}{Changepoint vector.} \item{k}{Number of changepoints.}
#' \item{network}{The network, a list of length NumSegs, where each element is
#' a NumNodes by NumNodes matrix.} \item{noise}{The standard deviation of the
#' applied Gaussian noise.}
#' @author Frank Dondelinger
#' @seealso \code{\link{generateNetwork}}
#' @examples
#' 
#' # Generate random network and simulate data with default parameters
#' dataset = simulateNetwork()
#' 
#' # Generate random network and simulate data with an average of 
#' # 1 change per node among network segments
#' dataset = simulateNetwork(lambda_3=1)
#' 
#' # Generate random network and simulate data with an average of 
#' # 1 change per node among network segments and standard deviation 
#' # of the Gaussian observation noise 0.5
#' dataset = simulateNetwork(lambda_3=1, noise=0.5)
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
#' @export simulateNetwork
simulateNetwork <-
function(l=100, min_phase_length=10, k_bar=10, q=10, 
  lambda_2=0.45, noise=0.25, net=NULL, lambda_3=2, spacing=0, 
  gauss_weights=FALSE, same=FALSE, changes='sequential', 
  fixed=FALSE, cps=NULL, saveFile=NULL) {

  if(is.null(net)) {
    net = generateNetwork(lambda_2, q, min_phase_length, k_bar, l, lambda_3,
                           spacing, gauss_weights, same, changes, fixed, cps)
  } else {
    l = net$l
    q = dim(net$network[[1]])[1]
  }
  
  network = net$network; epsilon = net$epsilon;
  k = net$k; changes = net$changes;
  
  # Simulate data from network
  sim_data = matrix(0, q, l) 
  sim_data[,1] = matrix(rnorm(q), q, 1);
  
  matrix_offset = 0; begin = 2;
  
  counter = 0;
  
  for (i in 1:(k+1)) {
	
	  parent_set = network[[i]];

	  change = epsilon[[i]];
	  
	  if(length(noise) == 1) {
	    seg_noise = noise;
	  } else {
	    seg_noise = noise[i];
	  }
	
    # Calculate regression model results for current segment  
    for (j in begin:change) {
	    new_pt = t(parent_set) %*% sim_data[,j-1]; 	
	  
	    # Nodes without parent are drawn from Gaussian
	    new_pt[new_pt == 0] = rnorm(sum(new_pt == 0));
	  
	    # Add noise
	    new_pt = new_pt + matrix(rnorm(q, 0, seg_noise), q, 1);
	  
	    # Scale to preserve variance = 1
	    sim_data[, j] = new_pt / sqrt(1 + seg_noise*seg_noise);
      
      # Fudge factor for correcting for instability  
      fudge = 10;
        
      if(any(abs(sim_data[,j]) > fudge*max(1, seg_noise))) {
        counter = counter + sum(abs(sim_data[,j]) > fudge*max(1, seg_noise));

        new_pt = sim_data[,j]; 
          
        new_pt[abs(new_pt) > fudge*max(1, seg_noise)] =
            new_pt[abs(new_pt) > fudge*max(1, seg_noise)] / fudge;
            
        sim_data[,j] = new_pt;
      }
	  }
	
	  begin = change + 1;
      
  }
  
  network.data = list(sim_data=sim_data, epsilon=epsilon, k=k, network=network, 
                 changes=changes, noise=noise)
  
  if(!is.null(saveFile))
    save(network.data, file=saveFile)
  
  return(network.data);
  
}


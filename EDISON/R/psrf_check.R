#' Check the potential scale reduction factors for all parameters (edges).
#' 
#' This function treats the edges of the network as parameters, calculates
#' their potential scale reduction factors and returns the highest value.
#' 
#' 
#' @param params Matrix of parameters.
#' @param q Number of nodes.
#' @param k_max Number of segments.
#' @param num_it Number of iterations/samples.
#' @return Returns the highest PSRF value.
#' @author Frank Dondelinger
#' @seealso \code{\link{psrf}}, \code{\link{psrf_check_hyper}}
#' @references Gelman and Rubin (1992) Inference from iterative simulation
#' using multiple sequences, Statistical Science.
#' @export psrf_check
psrf_check <-
function(params, q, k_max, num_it) {
  
  # Use only the last 5000 iterations
  limit = 5000;
  params_red = params[(num_it-limit):num_it,] != 0;
  
  PSRF_results = matrix(2, k_max, q)

  for(phase in 1:k_max) {
    for(gene in 1:q) {
      index = (q + 1)*(phase - 1) + gene
      param = params_red[,index];
      
      param_seqs = list()
      
      # Split sequence into 5 chains
      for(i in 1:5) {
        param_seqs[[i]] = t(matrix(param[((i-1)*(limit/5) + 1):(i*limit/5)]))
      }
      
      PSRF_results[phase, gene] = psrf(param_seqs)
      
    }
  }

  # Return best PSRF
  return(max(PSRF_results))
  
}


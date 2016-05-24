#' Checks the potential scale reduction factor for the hyperparameters.
#' 
#' This function checks the potential scale reduction factors for the
#' hyperparameters of the information sharing priors.
#' 
#' 
#' @param params Matrix of hyperparameters.
#' @param num_it Number of iterations/samples.
#' @return Returns the maximum PSRF value.
#' @author Frank Dondelinger
#' @seealso \code{\link{psrf}}, \code{\link{psrf_check}}
#' @references Gelman and Rubin (1992) Inference from iterative simulation
#' using multiple sequences, Statistical Science.
#' @export psrf_check_hyper
psrf_check_hyper <-
function(params, num_it) {
  # Only use last 5000 iterations
  limit = 5000;
  params_red = params[(num_it-limit):num_it,];
  
  PSRF_results = matrix(2, dim(params)[2], 1)

  for(p in 1:dim(params)[2]) {
    param = params_red[,p]
    
    param_seqs = list()
      
    # Split sequence into 5 chains
    for(i in 1:5) {
      param_seqs[[i]] = t(matrix(param[((i-1)*(limit/5) + 1):(i*limit/5)]))
    }
      
    PSRF_results[p] = psrf(param_seqs)
    
  }

  return(max(PSRF_results))
  
}


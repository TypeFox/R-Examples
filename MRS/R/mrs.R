#' Multi Resolution Scanning 
#'
#' This function executes the Multi Resolution Scanning algorithm to detect differences 
#' across multiple distributions.
#'
#' @param X Matrix of the data. Each row represents an observation.
#' @param G Numeric vector of the group label of each observation. Labels are integers starting from 1. 
#' @param n_groups Number of groups.
#' @param Omega Matrix defining the vertices of the sample space. 
#' The \code{"default"} option defines a hyperrectangle containing all the data points.
#' Otherwise the user can define a matrix  where each row represents a dimension,  
#' and the two columns contain the associated lower and upper limits for each dimension.
#' @param K Depth of the tree. Default is \code{K = 5}, while the maximum is \code{K = 14}.
#' @param init_state Initial state of the hidden Markov process.
#' The three states are \emph{null}, \emph{altenrative} and \emph{prune}, respectively.
#' @param beta Spatial clustering parameter of the transition probability matrix. Default is \code{beta = 1}.
#' @param gamma Parameter of the transition probability matrix. Default is \code{gamma = 0.3}.
#' @param eta Parameter of the transition probability matrix. Default is \code{eta = 0.3}.
#' @param alpha Pseudo-counts of the Beta random probability assignments. Default is \code{alpha = 0.5}.
#' @param return_global_null Boolean indicating whether to return the posterior probability of the global null hypothesis.
#' @param return_tree Boolean indicating whether to return the posterior representative tree. 
#' @param min_n_node Node in the tree is returned if there are more than \code{min_n_node} data-points in it.  
#' @return An \code{mrs} object. 
#' @export
#' @examples
#' set.seed(1) 
#' n = 20
#' p = 2
#' X = matrix(c(runif(p*n/2),rbeta(p*n/2, 1, 4)), nrow=n, byrow=TRUE)
#' G = c(rep(1,n/2), rep(2,n/2))
#' ans = mrs(X=X, G=G)
mrs <- function( X, 
                 G, 
                 n_groups = length(unique(G)), 
                 Omega = "default", 
                 K = 6, 
                 init_state = NULL, 
                 beta = 1.0, 
                 gamma = 0.3, 
                 eta = 0.3, 
                 alpha = 0.5,
                 return_global_null = TRUE,
                 return_tree = TRUE,
                 min_n_node = 0
                )
{
  X = as.matrix(X)
  if(Omega[1] == "default")
  {
    Omega = t(apply(X,2,range))
    Omega[,2] = Omega[,2]*1.0001
  }
  else
  {
    EmpOmega = t(apply(X,2,range))
    if( sum((EmpOmega[,1] >= Omega[,1]) + (EmpOmega[,2] < Omega[,2])) != 2*ncol(X)  )
    {
      print("ERROR: 'X' out of bounds, check 'Omega'")
      return(0);
    }
  }
  
  if( (K > 14) || (K < 0) )
  {
    print("ERROR: 0 <= K < 15")
    return(0);
  }
  
  if( min(G) < 1 )
  {
    print("ERROR: min(G) should be positive")
    return(0);
  }
  
  if(is.null(init_state)) {
    init_state = c((1-eta)*(1-gamma), (1-eta)*gamma, eta)
  } else {
    if( sum( init_state < 0 ) > 0 | sum( init_state ) != 1 )
    {
      print("ERROR: init_state should be non-negative and sum to 1")
      return(0);
    }
  }
  
  if( alpha<=0 )
  {
    print("ERROR: alpha > 0")
    return(0);
  }
  
    if( beta < 1 | beta > 2 )
    {
      print("ERROR: 1 <= beta <= 2")
      return(0);    
    }
  
  if( gamma < 0 | gamma > 1 )
  {
    print("ERROR: 0 <= gamma <= 1")
    return(0);
  }
  
  ans = fitMRScpp(X, G, n_groups, init_state, Omega, K, alpha, beta, gamma, eta, return_global_null, return_tree, min_n_node)

  if (return_tree) {
    ans$RepresentativeTree$EffectSizes = matrix( unlist(ans$RepresentativeTree$EffectSizes), 
                                                 nrow = length(ans$RepresentativeTree$Levels), byrow = TRUE)
    ans$RepresentativeTree$Regions = matrix( unlist(ans$RepresentativeTree$Regions), 
                                             nrow = length(ans$RepresentativeTree$Levels), byrow = TRUE)
  }

  
  colnames(ans$Data$X) = colnames(X)
  
  class(ans) = "mrs"
  return(ans)
}




SimulateDyadicLinearERGM <- function(N, dyadiccovmat, eta)
{
  
  # Note: Only performing an error check on the number of rows and columns
  # Not checking order of dyads or for duplication
  # For now this is probably OK; need to coordinate with function that builds the dyadic cov matrix
  
  invlogit <- function(y) return(exp(y)/(1+exp(y)))
  
  # do a little input checking
  n_par <- length(eta)
  n_dyad <- N*(N-1)/2
  if(dim(dyadiccovmat)[1] != n_dyad) stop("Invalid Input.")
  if (dim(dyadiccovmat)[2] != (n_par) + 2) 
    stop("Invalid Input. Need exactly one eta parameter for each non-id column of dyadic covariate matrix")
  
  # get log odds for each dyad
  if (length(eta)==1) dim(eta) <- c(1,1)
  logodds = dyadiccovmat[,2+(1:n_par)] %*% eta
  edge_prob = invlogit(logodds)
  
  # decide which dyads have edges present and return edge list matrix
  return(dyadiccovmat[edge_prob > runif(length(edge_prob)),1:2])
  
}

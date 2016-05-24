#' Calculates the potential scale reduction factor.
#' 
#' This function calculates the potential scale reduction factor of parameters
#' or hyperparameters over several MCMC simulations (or one simulation split
#' up). This can serve as a convergence diagnostic.
#' 
#' 
#' @param parameters A list of MCMC trajectories, where each trajectory is a
#' matrix with NumParams rows and NumIterations columns, where NumParams is the
#' number of parameters and NumIterations is the number of samples.
#' @return A vectors of length NumParams, containing the PSRF values for each
#' parameter.
#' @author Sophie Lebre
#' 
#' Frank Dondelinger
#' @seealso \code{\link{psrf_check}}, \code{\link{psrf_check_hyper}}
#' @references Gelman and Rubin (1992) Inference from iterative simulation
#' using multiple sequences, Statistical Science.
#' @examples
#' 
#' # Generate 5 'runs' of random samples from Gaussian N(0,1)
#' samples = list()
#' 
#' for(run in 1:5) {
#'   samples[[run]] = matrix(rnorm(1000), 1, 1000)
#' }
#' 
#' # Check potential scale reduction factor
#' # (Will be very close to 1 due to the samples being from 
#' # the same distribution)
#' psrf.val = psrf(samples)
#' 
#' 
#' # Now use slightly different Gaussian distributions for each 'run'.
#' for(run in 1:5) {
#'   mean = runif(1, 0, 2)
#'   samples[[run]] = matrix(rnorm(1000, mean, 1), 1, 1000)
#' }
#' 
#' # Check potential scale reduction factor
#' # (Should be > 1.1)
#' psrf.val = psrf(samples)
#' 
#' 
#' @export psrf
psrf <-
function(parameters) {
  
  # Number of sequences
  nbSeq=length(parameters)
  
  nbIterations=dim(parameters[[1]])[2]
  nbPars = dim(parameters[[1]])[1]
  
  # Compute B
  seq_means = matrix(0, nbPars, nbSeq)
  
  for(i in 1:nbSeq) {
    seq_means[, i] = apply(parameters[[i]], 1, mean)
  }
  
  B = nbIterations/(nbSeq-1)*
    apply((seq_means-matrix(
      apply(seq_means,1,mean),nbPars,nbSeq))^2,1,sum)
  
  # Compute W
  diffs = matrix(0, nbPars, nbSeq)
  
  for(i in 1:nbSeq) {
    diffs[, i] = apply((parameters[[i]] - kronecker(seq_means[, i], matrix(1, 1
                                                                           , nbIterations)))^2, 1, sum)
  }
  
  W = apply(diffs, 1, sum) / (nbSeq*(nbIterations-1))
  
  seq_overest = (nbSeq + 1) / nbSeq;
  
  PSRF = seq_overest*((nbIterations-1)/nbIterations+B/(W * nbIterations)) - 
    (nbIterations - 1)/(nbIterations*nbSeq)
  
  if(any(B == 0)) {
    PSRF[B==0] = seq_overest*((nbIterations-1)/nbIterations) - 
      (nbIterations - 1)/(nbIterations*nbSeq)
  }
    
  return(PSRF)
}


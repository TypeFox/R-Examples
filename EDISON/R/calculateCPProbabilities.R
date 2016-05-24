#' Calculate the changepoint probabilities.
#' 
#' This function calculates the marginal changepoint probabilities from the
#' changepoint samples taken during the MCMC simulation.
#' 
#' 
#' @param network.samples List of network and changepoint samples collected
#' during the MCMC simulation by \code{\link{EDISON.run}} and
#' \code{\link{runDBN}}.
#' @return Returns a matrix of dimension NumNodes by NumTimePoints, where each
#' entry contains the marginal posterior probability of a changepoint for that
#' node at that timepoint.
#' @author Frank Dondelinger
#' @examples
#' 
#' # Generate random gene network and simulate data from it
#' dataset = simulateNetwork()
#' 
#' # Run MCMC simulation to infer networks and changepoint locations
#' result = EDISON.run(dataset$sim_data, num.iter=500)
#' 
#' # Calculate posterior probabilities of changepoints
#' cps = calculateCPProbabilities(result)
#' 
#' @export calculateCPProbabilities
calculateCPProbabilities <-
function(network.samples) {
  timePoints = 1:(network.samples$n+1)
  
  sampled = network.samples[[1]]$sampled 
  numNodes = length(network.samples) - 1

  prob.cps = matrix(0, numNodes, length(timePoints)) 
  colnames(prob.cps) <- timePoints

  for(sample.i in 1:length(sampled)) {
     for(target in 1:numNodes) {
       cps.temp = network.samples[[target]]$cp_samples[sample.i,]
       cps.indices = which(timePoints %in% cps.temp)
       prob.cps[target,cps.indices] = prob.cps[target,cps.indices] + 1 
  
       numSegs = length(cps.temp) - 1
     } 
  }
  
  prob.cps = prob.cps / length(sampled) 
  
  global.cps = calculateCPPGlobal(prob.cps)
  
  return(list(node.cps=prob.cps, global.cps=global.cps))
}


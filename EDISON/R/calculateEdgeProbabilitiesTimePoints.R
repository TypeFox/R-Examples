#' Calculate the edge posterior probabilities for each timepoint.
#' 
#' This function calculates the marginal posterior edge probabilities of the
#' network at each timepoint.
#' 
#' 
#' @param network.samples Collection of network and changepoint samples of the
#' MCMC simulation, as obtained by \code{\link{EDISON.run}},
#' \code{\link{runDBN}}.
#' @param cps Changepoint vector.
#' @param numNodes Number of nodes in the network.
#' @return A list of length equal to the number of timepoints, where each entry
#' contains a matrix of size NumNodes by NumNodes with the marginal posterior
#' edge probabilities of the network at this timepoint.
#' @author Frank Dondelinger
#' @seealso \code{\link{calculateEdgeProbabilities}},
#' 
#' \code{\link{calculateEdgeProbabilitiesSegs}}
#' @export calculateEdgeProbabilitiesTimePoints
calculateEdgeProbabilitiesTimePoints <- 
  function(network.samples, cps, numNodes) {
  sampled = network.samples[[1]]$sampled 
  numSegs = length(cps) - 1
  
  segs = 2:cps[length(cps)]
  
  prob.networks = list() 
  
  for(seg in 1:length(segs)) {
    prob.networks[[seg]] = matrix(0, numNodes, numNodes)
  }
  
  for(sample.i in 1:length(sampled)) {
    for(target in 1:numNodes) {
      cps.temp = network.samples[[target]]$cp_samples[sample.i,]
      max.cp   = length(cps.temp) - 1
      cps.temp = cps.temp[cps.temp > 0]
      
      segs.temp = 2:cps.temp[length(cps.temp)]; seg = 1
      
      for(t in 1:length(segs.temp)) {
        
        if(t == cps.temp[seg+1]) {
          seg = seg + 1
        }
        
        segs.temp[t] = seg
      } 
      
      target.net.temp = 
        matrix(network.samples[[target]]$edge_samples[sample.i,],
               numNodes+1, max.cp)
      target.net.temp = (abs(target.net.temp) > 0)*1
      
      for(seg.i in 1:length(segs)) {
        prob.networks[[seg.i]][,target] = prob.networks[[seg.i]][,target] + 
          target.net.temp[1:numNodes, segs.temp[seg.i]] 
      }
    } 
  }
  
  for(seg.i in 1:length(segs)) {
    prob.networks[[seg.i]] = prob.networks[[seg.i]] / 
      length(sampled)
  }
  
  return(prob.networks)
}

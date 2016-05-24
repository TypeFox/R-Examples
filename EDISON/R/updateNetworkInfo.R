updateNetworkInfo <- function(network.info, Sall, Eall, prior.params=network.info$priorparams, 
                              posPhase=network.info$segment, target=network.info$target, 
                              q=network.info$q, self.loops=network.info$self.loops, 
                              k=network.info$k) {
  # Update network info datastructure (Not used in current implementation)
  
  if(any(Sall[[1]][Sall[[1]] != 0] != 1)) 
    stop('Illegal Network Structure in updateNetworkInfo.') 
  
  
  # If changepoint change, re-calculate global segments
  if(any(sapply(1:length(Eall), 
                function(target) 
                  Eall[[target]]!=network.info$cps[[target]]))) {
    network.info = CollectNetworkInfo(Sall, Eall, prior.params, posPhase, target, 
                                      q, self.loops, k)
  } 

  return(network.info)
}
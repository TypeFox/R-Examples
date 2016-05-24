#' Collects and saves output.
#' 
#' This function collects the network, changepoint and hyperparameter samples
#' taken from the MCMC simulation, and saves them to a file if appropriate.
#' 
#' 
#' @param counters List of counters for the number of moves that have been
#' proposed and accepted.
#' @param listStock Network, changepoint and hyperparameter samples.
#' @param GLOBvar Global variables of the MCMC simulation.
#' @param HYPERvar Hyperparameter variables.
#' @param OUTvar Output variables, including the output file.
#' @return Returns a list with an element for each target node which is also a
#' list. Each sublist containts the elements: \item{cp_samples }{Changepoint
#' samples, a NumSamples by MaxNumChangePoints matrix.}
#' \item{edge_samples}{Network samples (with regression parameters), a
#' NumSamples by (NumSegs * NumNodes) matrix.} \item{target}{The target node
#' for this subnetwork.} \item{hyper_samples}{Information sharing prior
#' hyperparameter samples, a NumSamples by NumHyperParams matrix.}
#' \item{sampled}{Sampled iterations.} \item{counters}{Counters for the number
#' of proposed and accepted moves.}
#' @author Frank Dondelinger
#' @export output
output <-
function(counters, listStock, GLOBvar, HYPERvar, OUTvar){

  ### Assignment of global variables used here ###
  target = GLOBvar$target
  n = GLOBvar$n
  q = GLOBvar$q
  ### End assignement ###

  ### Assignment of results variables used here ###
  # Counters
  cptMove = counters$cptMove
  acceptMove = counters$acceptMove
  # ListStock
  Estock = listStock$Estock
  Bstock = listStock$Bstock
  hyperstock = listStock$hyperstock
  ### End assignment ###

  ### Assignment of output variables used here ###
  outputFile=OUTvar$outputFile
  analysis = OUTvar$analysis
  ### End assignment ###
 
  results.all = list()

  n_samples = dim(Bstock[[1]])[1]
  after.burnin = round(n_samples/4):n_samples

  # Sample starting from after burnin (default: 1/4 of run)
  if(length(after.burnin) > 1000) {
    sampled = sort(sample(after.burnin, 1000))
  } else {
    sampled = after.burnin
  }

  for(target in 1:q) {
              
    # Collect results
    results = list(cp_samples=Estock[[target]][sampled,], 
      edge_samples = Bstock[[target]][sampled,], 
      target=target, hyper_samples=hyperstock[sampled,], 
      sampled=sampled, counters=counters)
        
    results.all[[target]] = results 
  
    if(OUTvar$by.node && OUTvar$save.file) {
      save(results, file=paste(outputFile, "_analysis_", target, 
	                       sep=""))
    }
  }
      
  if(!OUTvar$by.node && OUTvar$save.file) {
    results = results.all
    save(results, file=paste(outputFile, "_analysis", 
	     sep=""))
	}
	
  results.all$n = n

  return(results.all)
}


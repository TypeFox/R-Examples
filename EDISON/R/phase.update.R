#' Make a network structure or hyperparameter move.
#' 
#' This function makes a network structure or information sharing
#' hyperparameter move.
#' 
#' 
#' @param Eall List of changepoints with one entry for each target node. Each
#' entry has length equal to the number of changepoints for that target node.
#' @param Sall Network structure: List of length equal to the number of target
#' nodes, where each list entry is a NumSegs by NumNodes matrix.
#' @param Ball Network structure with regression coefficients: Same as Sall,
#' but with regression coefficients as matrix entries.
#' @param Sig2all Sigma squared.
#' @param X Input response data.
#' @param Y Input target data.
#' @param GLOBvar Global variables used during the MCMC simulation.
#' @param HYPERvar Hyperparameter variables.
#' @param target Current target node.
#' @return Returns a list with the following elements: \item{E}{Changepoints
#' for the current target node.} \item{Sall}{Network structure (possibly
#' updated).} \item{Ball}{Network structure regression coefficients (possibly
#' updated).} \item{Sig2all}{Sigma squared.} \item{prior.params}{Information
#' sharing prior hyperparameters (possibly updated).} \item{k}{Level-2
#' exponential prior hyperparameter (possibly updated).} \item{move}{Move type:
#' 4 for a network structure move, 5 hyperparameter move.}
#' \item{move}{Structure Move type: 1 for a network structure move, 2 for a
#' level-1 hyperparameter move, 3 for a level-2 hyperparameter move.}
#' \item{accept}{1 if the move has been accepted, 0 otherwise.}
#' @author Sophie Lebre
#' 
#' Frank Dondelinger
#' @seealso \code{\link{make_structure_move}}
#' @references For more information on network structure moves and information
#' sharing priors, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export phase.update
phase.update <-
function(Eall, Sall, Ball, Sig2all, X, Y, GLOBvar, HYPERvar, 
  target) {
  
  E = Eall[[target]]
  # Current number of changepoints
  s = length(E) - 2
  
  ### Assignment of global variables used here ###
  q = GLOBvar$q
  qmax = GLOBvar$qmax
  Mphase = GLOBvar$Mphase
  nbVarMax = GLOBvar$nbVarMax
  smax = GLOBvar$smax
  lmax = GLOBvar$lmax
  method = GLOBvar$method
  small_prop = GLOBvar$small_prop
  self.loops = GLOBvar$self.loops
  fixed.edges = GLOBvar$fixed.edges
  ### End assignment ###

  ### Assignment of hyperparameters variables used here ###
  c = HYPERvar$c
  alphalbd = HYPERvar$alphalbd
  betalbd = HYPERvar$betalbd
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  k = HYPERvar$k
  ### End assignment ###

  prior.params = HYPERvar$prior.params
        
  Sig2 = Sig2all[[target]]
    
  ## Observations in the chosen phase
  y = Y[[target]]
  x = X[[target]]

  model = 0
    
  # Group information about network segments    
  network.info = CollectNetworkInfo(Sall, Eall, prior.params, -1, 
                                    target, q, self.loops, k)
    
  # Try hyperparameter move
  hyper.move = HyperparameterMove(method, network.info, GLOBvar, 
                                  HYPERvar$hyper.proposals)
    
  # If no move made, try structure move
  if(!hyper.move$move.made) {
    ## Compute a structure move 
    bduout =  make_structure_move(x, y, Sall[[target]], Ball[[target]], Sig2, q, 
                qmax, network.info, method, Mphase, E, fixed.edges[,target],
                HYPERvar)
			
    Sall[[target]] = bduout$newS
    Ball[[target]] = bduout$newB

  } else {
    # Hyperparameter move
    prior.params = hyper.move$network.info$prior.params
    k = hyper.move$network.info$k
    bduout = list(move=hyper.move$move, accept=hyper.move$accept)
    
  }
  
  # Update the regression coefficient
  for(posPhase in 1:(s+1)) {  
    ## Update Sig2: case MultiVar
    if(nbVarMax >1){
      Sig2 = Sig2all[[target]][posPhase]
      Sig2all[[target]][posPhase] = 
        updateSigMulti(E[posPhase], X[[target]], Y[[target]], E, 
          Sall[[target]], Ball[[target]], Sig2, Mphase, alphad2, betad2, 
                                                   v0, gamma0)
    } #end update Sig2
  } # end update each phase

  ## Update Sig2: case UniVar
  if(nbVarMax == 1){
    for(target in 1:q) {
      Sig2all[[target]] = 
        updateSigSolo(X[[target]], Y[[target]], E, Sall[[target]], 
                            Ball[[target]], Sig2all[[target]], Mphase, alphad2, 
                            betad2, v0, gamma0)
    }
  }
  
  ##  Return all variables
  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return( list( E=E, Sall=Sall, Ball=Ball, Sig2all=Sig2all, 
                prior.params=prior.params, k=k,
                accept=bduout$accept, move=bduout$move))
}


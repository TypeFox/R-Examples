#' Make changepoint death move.
#' 
#' This function makes a changepoint death move, possibly removing a
#' changepoint.
#' 
#' 
#' @param Eall Changepoints: List of target nodes, where each element contains
#' a vector of changepoints.
#' @param Sall Network structure: List of target nodes, where each element is a
#' NumSegs by NumNodes matrix giving the parents for the target node in each
#' segment. A binary matrix.
#' @param Ball Network parameters: Similar to network structure, but with
#' regression parameters included.
#' @param Sig2all Sigma squared parameters.
#' @param X Response data.
#' @param Y Target data.
#' @param D Hyperparameter.
#' @param GLOBvar Global variables of the MCMC simulation.
#' @param HYPERvar Hyperparameter variables.
#' @param target Which target node the move is being proposed for.
#' @return A list with elements: \item{E}{New changepoint vector for target
#' node.} \item{Sall}{Updated network structure.} \item{Ball}{Updated network
#' structure with regression parameters.} \item{Sig2all}{Updated sigma
#' squared.} \item{prior.params}{Updated vector of structure prior
#' hyperparameters.} \item{accept}{Whether the move was accepted or not.}
#' \item{move}{What type of move was made. In this case \code{move=2} for a
#' changepoint death move.} \item{alpha}{The acceptance ratio of the move.}
#' \item{estar}{The location of the removed changepoint.}
#' \item{k}{Hyperparameter.}
#' @author Sophie Lebre
#' 
#' Frank Dondelinger
#' @seealso \code{\link{cp.birth}}, \code{\link{cp.shift}}
#' @references For more information about the different changepoint moves, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export cp.death
cp.death <-
function(Eall, Sall, Ball, Sig2all, X, Y, D, GLOBvar, HYPERvar,
                     target){
  ### INPUT:  E, Sall, Ball, Sig2all, X, Y, D, GLOBvar, HYPERvar
  ### OUTPUT: 

  E = Eall[[target]]

  # Current number of changepoints
  s = length(E) - 2
  
  ### Assignment of global variables used here ###
  q = GLOBvar$q
  Mphase = GLOBvar$Mphase
  nbVarMax = GLOBvar$nbVarMax
  smax = GLOBvar$smax
  lmax = GLOBvar$lmax
  method = GLOBvar$method
  ## End assignment ###

  ### Assignment of hyperparameters variables used here ### 
  alphalbd = HYPERvar$alphalbd
  betalbd = HYPERvar$betalbd
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  prior.params = HYPERvar$prior.params;
  k.par = HYPERvar$k
  ### End assignment ###
   
 
  ## Sample the CP to be removed
  estar = sample(c(E[2:(length(E)-1)], E[2:(length(E)-1)]), 1)
  
  ##  Position of the phase starting at the selected CP
  poskstar = sum(E <= estar)
  
  # Sample the edge vector to be maintained  for merge phase after 
  # removing the CP
  # (newRight =1 if the edge vector to the Right of the CP is maintained, 
  # =0 otherwise)
  newRight = sample(0:1,1)
  
  ## Position of the line of the phase to be removed in the matrices Sall and Ball
  away = (1-newRight) * poskstar + newRight * (poskstar-1)
  
  # Compute the matrices required for the computation of the 
  # acceptance probability alpha
  yL = Y[(Mphase[E[poskstar-1]]:(Mphase[estar]-1))]
  xL = X[(Mphase[E[poskstar-1]]:(Mphase[estar]-1)),]
  yR = Y[(Mphase[estar]:(Mphase[E[poskstar+1]]-1))]
  xR = X[(Mphase[estar]:(Mphase[E[poskstar+1]]-1)),]
  y2 = array(c(yL,yR))
  x2 = rbind(xL,xR)
 
  if(nbVarMax>1){
    Sig2 = Sig2all[poskstar - 1 + newRight]
  } else {
    Sig2 = Sig2all
  }

  ## Update delta
  delta2 = sampleDelta2(poskstar-1+newRight, x2, q, 
                        Ball[[target]], Sall[[target]], Sig2, 
                        alphad2, betad2)

  # Compute projection of the matrices required for the computation 
  # of the acceptance probability alpha
  Px2 = computePx(length(y2), 
                  as.matrix(x2[,
                    which(Sall[[target]][poskstar-1+newRight,] == 1)]), 
                  delta2) 
  PxL = computePx(length(yL), 
                  as.matrix(xL[,
                    which(Sall[[target]][poskstar-1,] == 1)]), 
                  delta2)
  PxR = computePx(length(yR), 
                  as.matrix(xR[,
                    which(Sall[[target]][poskstar,] == 1)]), 
                  delta2)

  prior_ratio = 1;
  proposal.ratio = 1;

  remove_right = !newRight
    
  if(method != 'poisson') {
    
    # Proposal Ratio Under Poisson Proposals
    s.away = sum(Sall[[target]][away,1:q])
    lambda = rgamma(1, shape=s.away + alphalbd, rate=1 + betalbd) 
    
    p.old = (factorial(q-s.away)/factorial(q)) * lambda ^ s.away
    proposal.ratio = p.old

    # Prior ratio
    network.info.old = 
      CollectNetworkInfo(Sall, Eall, prior.params, -1, target, q,
                         -1, k.par)
   
    Sall.new = Sall 
    Sall.new[[target]] = 
      matrix(Sall[[target]][(1:(s+1))[-c(away)],], s, q+1) 
    E.new = E[(1:(s+2))[-c(poskstar)]]
    
    Eall.new = Eall
    Eall.new[[target]] = E.new

    network.info.new = CollectNetworkInfo(Sall.new, Eall.new, prior.params, -1, 
                         target, q,
                         -1, k.par)

    if(method == 'exp_hard' || method == 'exp_soft') {
      log.prob.old = NetworkProbExp(network.info.old)
      log.prob.new = NetworkProbExp(network.info.new)
    } else if(method == 'bino_hard' || method == 'bino_soft') {
      log.prob.old = NetworkProbBino(network.info.old)
      log.prob.new = NetworkProbBino(network.info.new)
    } 
    
    prior_ratio = exp(log.prob.new - log.prob.old)
   
  }
 
  pp.ratio = (prior_ratio*proposal.ratio)

  ## Compute the acceptance probability alpha
  alpha = bp.computeAlpha(-1, sum(Sall[[target]][poskstar-1+newRight,])-1, s-1, 
                   Mphase[E[poskstar-1]], Mphase[estar], Mphase[E[poskstar+1]], 
                   yL, PxL, yR, PxR, y2, Px2, D, delta2, q, smax, v0, gamma0, 
                   1/pp.ratio)
  
  ## Sample u to decide on acceptance
  u = runif(1,0,1)

  # Boolean for the acceptance of the CP death move initially set to 0 
  # (=1 if birth accepted, 0 otherwise)
  accept = 0
  
  if(!is.nan(alpha) & u <= alpha){
    ## Acceptance of the death of the selected CP
    
    accept=1
   
    ## Remove the CP in E and the phase in the matrices Sall and Ball
    if(nbVarMax>1){
      Sig2all = Sig2all[(1:(s+1))[-c(away)]]
    }
    
    E = E[(1:(s+2))[-c(poskstar)]]
    Sall[[target]] = matrix(Sall[[target]][(1:(s+1))[-c(away)],], s, q+1)
    Ball[[target]] = matrix(Ball[[target]][(1:(s+1))[-c(away)],], s, q+1)
  }

  #  Return all variables
  # (+ variable move describing the move type  
  # (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return( list( E=E, Sall=Sall, Ball=Ball, Sig2all=Sig2all, 
                prior.params=prior.params, accept=accept, move=2, 
                alpha=alpha, estar=estar, k=k.par))
}


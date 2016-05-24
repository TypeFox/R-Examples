#' Make changepoint birth move.
#' 
#' This function makes a changepoint birth move, possibly adding a changepoint.
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
#' structure with regression parameters.} \item{Sig2all}{Udated sigma squared.}
#' \item{prior.params}{Updated vector of structure prior hyperparameters.}
#' \item{accept}{Whether the move was accepted or not.} \item{move}{What type
#' of move was made. In this case \code{move=1} for a changepoint birth move.}
#' \item{alpha}{The acceptance ratio of the move.} \item{estar}{The location of
#' the new changepoint.} \item{k}{Hyperparameter.}
#' @author Sophie Lebre
#' 
#' Frank Dondelinger
#' @seealso \code{\link{cp.death}}, \code{\link{cp.shift}}
#' @references For more information about the different changepoint moves, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export cp.birth
cp.birth <-
function(Eall, Sall, Ball, Sig2all, X, Y, D, 
                     GLOBvar, HYPERvar, target){
  # INPUT: E, Sall, Ball ,Sig2all, X, Y, D, GLOBvar, HYPPERvar
  # OUTPUT: Result of move
  
  E = Eall[[target]]
  
  # Current number of changepoints
  s = length(E) - 2
  
  ### Assignment of global variables used here ###
  q = GLOBvar$q
  qmax = GLOBvar$qmax
  Mphase = GLOBvar$Mphase
  minPhase = GLOBvar$minPhase
  nbVarMax = GLOBvar$nbVarMax
  smax = GLOBvar$smax
  dyn = GLOBvar$dyn
  lmax = GLOBvar$lmax
  small_prop = GLOBvar$small_prop
  method = GLOBvar$method 
  self.loops = GLOBvar$self.loops
  fixed.edges = GLOBvar$fixed.edges[,target]
  ### End assignment ###

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
 
  S = Sall[[target]]
  B = Ball[[target]]
 
  # Search for possible CP, not in E and not close to E if 
  # minPhase (length of phase) is > than 1
  toremove = E
  
  if(minPhase>1) 
    for(i in 1:(minPhase-1)) 
      toremove = c(toremove, E-i, E+i)
  
  # Possible CPs are those not in 'toremove'
  possibleCP = setdiff((1+dyn):E[length(E)], toremove)
  
  # Sample the new CP "estar"
  estar = sample(c(possibleCP, possibleCP),1)
  E.new = sort(c(E, estar))

  # Position of the phase containing the new CP
  poskl = sum(E < estar)
  
  # Current edges vector S in the phase containing the new CP
  Sold = S[poskl,]
  
  # Current number of edges k in the phase containing the new CP
  k = sum(Sold) - 1

  # Sample lambda
  lambda = rgamma(1, shape=alphalbd, rate=betalbd)

  # Sample a new edges vector newS
  newS = array(1, q+1)

  # Sample Right or Left (the position for the new 
  # edges vector: to the right or to the left of the new CP estar)
  newRight = runif(1,0,1) >= 1/2
	
  newS[1:q] = 1:q %in% sample(1:q, sampleK(0, qmax, lambda, 1), replace=FALSE)
  
  if(!newRight){
    # New edges vector to the left of the new CP
    # Boolean (= 1 if  the new edges vector is to the right 
    # of the new CP, 0 otherwise) 
    sL = newS
    sR = Sold
  } else {
    # New edges vector to the right of the new CP
    # Boolean (= 1 if  the new edges vector is to the right of 
    # the new CP, 0 otherwise)
    sR = newS
    sL = Sold
  }
  
  # Compute the matrices required for the computation of the 
  # acceptation probability alpha
  yL = Y[(Mphase[E[poskl]]:(Mphase[estar]-1))]
  xL = X[(Mphase[E[poskl]]:(Mphase[estar]-1)),]
  yR = Y[(Mphase[estar]:(Mphase[E[poskl+1]]-1))]
  xR = X[(Mphase[estar]:(Mphase[E[poskl+1]]-1)),]
  y2 = array(c(yL, yR))
  x2 = rbind(xL, xR)
  
  ## Updating parameters
  if(nbVarMax > 1){
    Sig2 = Sig2all[poskl]
  } else {
    Sig2 = Sig2all
  }

  delta2 = sampleDelta2(poskl, x2, q, B, S, Sig2, alphad2, betad2)
  
  # Compute projection of the matrices required for the computation 
  # of the acceptation probability alpha
  PxL = computePx(length(yL), as.matrix(xL[,which(sL == 1)]), delta2)
  PxR = computePx(length(yR), as.matrix(xR[,which(sR == 1)]), delta2)
  Px2 = computePx(length(y2), as.matrix(x2[,which(Sold == 1)]), delta2)
  
  prior_ratio = 1;
  proposal.ratio = 1; 
  
  # Calculate information sharing priors (if applicable)
  if(method != 'poisson') {
    # Proposal Ratio Under Poisson Proposals
    s.new = sum(newS[1:q])  

    p.new = (factorial(q-s.new)/factorial(q)) * lambda ^ s.new
    proposal.ratio = 1 / p.new

    # Prior Ratio
    network.info.old = 
      CollectNetworkInfo(Sall, Eall, prior.params, -1, target, q,
                         -1, k.par)
    
    # Create proposed new network
    Sall.new = Sall
    Sall.new[[target]] = 
      matrix(0, dim(Sall.new[[target]])[1] + 1, 
                                   dim(Sall.new[[target]])[2])
    Sall.new[[target]][1:poskl,] = Sall[[target]][1:poskl,]
    
    Sall.new[[target]][(poskl + 1):(s+2),] = Sall[[target]][(poskl:(s+1)),]      
      
    if(newRight) {
      Sall.new[[target]][poskl+1,] = newS
    } else {
      Sall.new[[target]][poskl,] = newS
    }
    
    Eall.new = Eall
    Eall.new[[target]] = E.new

    network.info.new = 
      CollectNetworkInfo(Sall.new, Eall.new, prior.params, -1, 
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

  
  ## Compute the acceptance probability alpha
    
  pp.ratio = (prior_ratio*proposal.ratio) 

  alpha = bp.computeAlpha(1, sum(newS)-1, s, Mphase[E[poskl]], 
                                  Mphase[estar], Mphase[E[poskl+1]], yL, 
                                  PxL, yR, PxR, y2, Px2, D, delta2, q, 
                                  smax, v0, gamma0, pp.ratio)
  			 
  ## Sample u to conclude either to  acceptation or to rejection
  u = runif(1,0,1)
  
  # Boolean for the acceptance of the CP birth move initially set 
  # to 0 (=1 if birth accepted, 0 otherwise)
  accept = 0

  if(!is.nan(alpha) & u <= alpha && 
    AcceptableMove(t(as.matrix(newS)), qmax, self.loops, 
                                      target, fixed.edges)) {
   
    ## Acceptance of the birth of the new CP
    ## Move acceptance boolean =1
    accept = 1

    ## Compute new Sig2
    if(nbVarMax > 1){
      newSig2 = array(0, s+2)
      newSig2[(1:(s+2))[-c(poskl,poskl+1)]] = Sig2all[(1:(s+1))[-c(poskl)]]
    }
    
    ## Compute new regression parameters newB
    newB = matrix(0, s+2, q+1)
    newB[(1:(s+2))[-c(poskl,poskl+1)],] = B[(1:(s+1))[-c(poskl)],]  
	
    ## Update newB 
    if(!newRight){
      ## Update the phase to the left of the new CP (in newB)
      if(nbVarMax > 1){
        newSig2[poskl] = sampleSig2(yL,PxL,v0,gamma0)
        newSig2[poskl+1] = Sig2all[poskl]
        Sig2all = newSig2
        Sig2 = newSig2[poskl]
      }
      
      newB[poskl+1,] = B[poskl,]
      newB[poskl, which(newS == 1)] = 
        sampleBxy(xL[,which(newS==1)], yL, Sig2, delta2)
	  
    } else {
      ## Update the phase to the right of the new CP (in newB)
      if(nbVarMax > 1){
        newSig2[poskl+1] = sampleSig2(yR, PxR, v0, gamma0)
        newSig2[poskl] = Sig2all[poskl]
        Sig2all = newSig2
        Sig2 = newSig2[poskl]
      }
      
      newB[poskl,] = B[poskl,]
      newB[poskl+1, which(newS == 1)] = 
        sampleBxy(xR[, which(newS == 1)], yR, Sig2, delta2)
	  
    }
    
    ## Update current model and parameters
    Ball[[target]] = newB
    Sall[[target]] = (abs(Ball[[target]])>0)*1
    Eall[[target]] = E.new
  }

  #  Return all variables
  # (+ variable move describing the move type  
  # (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return(list(E=Eall[[target]], Sall=Sall, Ball=Ball, Sig2all=Sig2all, 
              prior.params=prior.params, accept=accept, move=1, 
              alpha=alpha, estar=estar, k=k.par))
}


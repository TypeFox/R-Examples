#' Makes a changepoint shift move.
#' 
#' This function makes a changepoint shift move, possibly moving one of the
#' changepoints.
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
#' @seealso \code{\link{cp.birth}}, \code{\link{cp.death}}
#' @references For more information about the different changepoint moves, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export cp.shift
cp.shift <-
function(Eall, Sall, Ball, Sig2all, X, Y, GLOBvar, 
                     HYPERvar, target){
  ### INPUT: u,rho,E,S,B,Sig2,X,Y
  ### OUTPUT: 
  
  E = Eall[[target]]
  
  delta2 = 1
  
  # Current number of changepoints
  s = length(E) - 2
  
  ### Assignment of global variables used here ###
  q = GLOBvar$q
  Mphase = GLOBvar$Mphase
  minPhase = GLOBvar$minPhase
  nbVarMax = GLOBvar$nbVarMax
  smax = GLOBvar$smax
  ## End assignment ###

  ### Assignment of hyperparameters variables used here ###  
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  betas = HYPERvar$betas
  k.par = HYPERvar$k
  ### End assignment ###
  
  prior.params = HYPERvar$prior.params
  
  ## Sample the CP to be shifted
  estar = sample(c(E[2:(length(E)-1)], E[2:(length(E)-1)]), 1)

  ## Position of the phase starting at the selected CP
  poskstar = sum(E <= estar)

  ## Possible new position for the selected CP (CP-1, CP+1)
  newCPs = c(E[poskstar]-1,E[poskstar]+1)
  
  # Remove positions that create too short phases (minPhase)
  newCPs = newCPs[which( !( newCPs %in% c(E[poskstar-1],
                                          E[poskstar-1]+minPhase-1, 
                                          E[poskstar+1],
                                          E[poskstar+1]-minPhase+1,E)))]
  
  # Boolean for the acceptation of the CP shift move initially 
  # set to 0 (=1 if birth accepted, 0 otherwise)
  accept = 0
  
  ## If there is at least one option to shift the selected CP 
  if(length(newCPs) > 0){
    ## Sample new CP position
    newCP = sample( c(newCPs, newCPs), 1)
    
    # Compute the matrices required for the computation of the 
    # acceptance probability alpha
    
    # Matrices for the model before shifting the CP
    yL = Y[ Mphase[E[poskstar-1]]:(Mphase[estar]-1) ]
    xL = X[ Mphase[E[poskstar-1]]:(Mphase[estar]-1), ]
    yR = Y[ Mphase[estar]:(Mphase[E[poskstar+1]]-1) ]
    xR = X[ Mphase[estar]:(Mphase[E[poskstar+1]]-1), ]

	
    ## Matrices for the model after shifting the CP
    yLStar = Y[ Mphase[E[poskstar-1]]:(Mphase[newCP]-1) ]
    xLStar = X[ Mphase[E[poskstar-1]]:(Mphase[newCP]-1), ]
    yRStar = Y[ Mphase[newCP]:(Mphase[E[poskstar+1]]-1) ]
    xRStar = X[ Mphase[newCP]:(Mphase[E[poskstar+1]]-1), ]
    
    if(nbVarMax > 1){
      
      ## Update delta for the left phase
      delta2 = sampleDelta2(poskstar-1, xL, q, Ball, Sall, 
                            Sig2all[poskstar-1], alphad2, betad2)
      
      ## Compute projection of the matrices for the left phase
      PxL = computePx(length(yL), as.matrix(xL[, which(Sall[poskstar-1,] == 1)]), 
                      delta2)
      PxLStar = computePx(length(yLStar), 
                          as.matrix(xLStar[, which(Sall[poskstar-1,] == 1)]), 
                          delta2)
      
      ## Update delta for the right phase     
      delta2 = sampleDelta2(poskstar, xR, q, Ball, Sall, Sig2all[poskstar], 
                            alphad2, betad2)
      
      ## Compute projection of the matrices for the right phase
      PxR = computePx(length(yR), as.matrix(xR[, which(Sall[poskstar,] == 1)]), 
                      delta2)
      PxRStar = computePx(length(yRStar), 
                          as.matrix(xRStar[, which(Sall[poskstar,] == 1)]), 
                          delta2)
    } else {
      ## Update delta for the left phase
      delta2 = sampleDelta2(poskstar-1, xL, q, Ball, Sall, Sig2all, 
                            alphad2, betad2)
      ## Compute projection of the matrices for the left phase
      PxL = computePx(length(yL), as.matrix(xL[, which(Sall[poskstar-1,] == 1)]), 
                      delta2)
      PxLStar = computePx(length(yLStar), 
                          as.matrix(xLStar[, which(Sall[poskstar-1,] == 1)]), 
                          delta2)
      
      ## Update delta for the right phase     
      delta2 = sampleDelta2(poskstar, xR, q, Ball, Sall, Sig2all, 
                            alphad2, betad2)
      
      ## Compute projection of the matrices for the right phase
      PxR = computePx(length(yR), 
                      as.matrix(xR[, which(Sall[poskstar,] == 1)]), 
                      delta2)
      PxRStar = computePx(length(yRStar), 
                          as.matrix(xRStar[, which(Sall[poskstar,] == 1)]), 
                          delta2)
    }

    ## Compute the logarithm of the Likelihood Ratio (LR)
    logLR=lgamma(((Mphase[newCP]-Mphase[E[poskstar-1]])+v0)/2)+
      lgamma(((Mphase[E[poskstar+1]]-Mphase[newCP])+v0)/2)-
      lgamma(((Mphase[E[poskstar]]-Mphase[E[poskstar-1]])+v0)/2)- 
      lgamma(((Mphase[E[poskstar+1]]-Mphase[E[poskstar]])+v0)/2)+ 
      ((Mphase[E[poskstar]]-Mphase[E[poskstar-1]])+v0)/2* 
      log((gamma0+t(yL)%*%PxL%*%yL)/2)+
      ((Mphase[E[poskstar+1]]-Mphase[E[poskstar]])+v0)/2*
      log((gamma0+t(yR)%*%PxR%*% yR)/2)-
      (((Mphase[newCP]-Mphase[E[poskstar-1]])+v0)/2)*
      log((gamma0+t(yLStar)%*%PxLStar %*%yLStar)/2)-
      ((Mphase[E[poskstar+1]]-Mphase[newCP])+v0)/2*
      log((gamma0+t(yRStar) %*% PxRStar%*%yRStar)/2)
  

    ## New CP vector Estar
    Estar = E
    Estar[poskstar] = newCP

    # Computation of the proposal Ratio
    # Vector of length the current number of phases = 
    # c(1,2,2,...,2,2,1) i.e. the number of CP that can 
    # potentially be shifted into each phase
    nbmove = c(1,array(2,s-1),1)
    propRatio = (2*s-sum(((E[2:(s+2)]-E[1:(s+1)]) <= minPhase) * nbmove))/
      (2 * s - sum(((Estar[2:(s+2)]-Estar[1:(s+1)]) <= minPhase) * nbmove))

    ## Computation of alpha
    if(!is.nan(logLR) & (logLR+log(propRatio))<0){
      alpha = min(c(1, exp(logLR)*propRatio))
    } else {
      alpha = 1
    }

    ## Sample u to decide whether the CP shift is accepted or not
    u = runif(1,0,1)

    # Boolean for the acceptance of the CP move
    # (=1 if birth accepted, 0 otherwise)

    if(u <= alpha){
      ## Acceptance of the move of the selected CP
      ## Move acceptation boolean =1
      accept = 1
      E[poskstar] = newCP
    }
  }
  
  Eall[[target]] = E
  
  # Return all variables
  # (+ variable move describing the move type  
  # (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return(list(E=E, Sall=Sall, Ball=Ball, Sig2all=Sig2all, 
              prior.params=prior.params, accept=accept, move=3, 
              k=k.par))
}


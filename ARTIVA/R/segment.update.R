segment.update <-
function(E, Sall, Ball, Sig2all, X, Y, GLOBvar, HYPERvar, cptMove, acceptMove){
  
  # current number of changepoints
  s = length(E) - 2
  
  ### assignement of global variables used here ###
  q = GLOBvar$q
  qmax = GLOBvar$qmax
  Mphase = GLOBvar$Mphase
  nbVarMax = GLOBvar$nbVarMax
  smax = GLOBvar$smax
  ### end assignement ###

  ### assignement of hyperparameters variables used here ###
  c = HYPERvar$c
  alphalbd = HYPERvar$alphalbd
  betalbd = HYPERvar$betalbd
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  ### end assignement ###
 
  ## For each current phase
  for (phase in E[1:(s+1)]){
  
    ## Position of the phase
    posPhase = which(E==phase)

    ## Parameters
    B = Ball[posPhase,]
    S = (abs(B) > 0) * 1
    k = sum(S)-1
    if(nbVarMax >1){
      Sig2 = Sig2all[posPhase]
    } else {
      Sig2 = Sig2all
    }
    
    ## Observations in the chosen phase
    y = Y[ Mphase[phase]:(Mphase[E[posPhase+1]]-1) ]
    x = X[ Mphase[phase]:(Mphase[E[posPhase+1]]-1), ]
    
    ## Updating hyperparameters
    delta2 = rinvgamma(1, shape=k + alphad2, scale=betad2 + B[which(S==1)] %*% t(x[,which(S==1)]) %*% x[,which(S==1)] %*% B[which(S==1)] / (2*Sig2) ) 
    lambda = rgamma(1, shape=k + alphalbd, rate=1 + betalbd)

    ## Compute acceptation probability vector rho
    rho3 = computeRho3(k, 0, qmax, c, lambda)
   
    ## Sample u
    u = runif(1, 0, 1)

    ## Compute the corresponding move (Edge birth, Edge death or Update the regression coefficient) 
    bduout = bdu(u, rho3, x, y, S, Sig2, delta2, q, v0, gamma0)

    cptMove[bduout$move]=cptMove[bduout$move] + 1
    acceptMove[bduout$move]=acceptMove[bduout$move] + bduout$accept
    
    Sall[posPhase,] = bduout$newS
    Ball[posPhase,] = bduout$newB

    ## Update Sig2
    if(nbVarMax >1){
      Sig2all[posPhase] = updateSigMulti(phase, X, Y, E, Sall, Ball, Sig2, Mphase, alphad2, betad2, v0, gamma0)
    } #end update Sig2
  } # end update each phase

  ##  Return all variables
  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return(list(E=E, Sall=Sall, Ball=Ball, Sig2all=Sig2all, accept=1, move=4, EdgesMove=cptMove, EdgesAccept=acceptMove))

}

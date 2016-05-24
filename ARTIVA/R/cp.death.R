cp.death <-
function(E, Sall, Ball, Sig2all, X, Y, D, GLOBvar, HYPERvar){
  ### INPUT:  E, Sall, Ball, Sig2all, X, Y, D, GLOBvar, HYPERvar
  ### OUTPUT: 
  ### depends on: .

  # current number of changepoints
  s = length(E) - 2
  
  ### assignement of global variables used here ###
    m=GLOBvar$m
	q = GLOBvar$q
  Mphase = GLOBvar$Mphase
  nbVarMax = GLOBvar$nbVarMax
  smax = GLOBvar$smax
  ### end assignement ###

  ### assignement of hyperparameters variables used here ###
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  ### end assignement ###
  
  ## Sample the CP to be removed
  estar = sample(c(E[2:(length(E)-1)], E[2:(length(E)-1)]), 1)
  ##  Position of the phase starting at the selected CP
  poskstar = sum(E <= estar)

  ## Compute the matrices required for the computation of the acceptation probability alpha
  yL = Y[(Mphase[E[poskstar-1]]:(Mphase[estar]-1))]
  xL = X[(Mphase[E[poskstar-1]]:(Mphase[estar]-1)),]
  yR = Y[(Mphase[estar]:(Mphase[E[poskstar+1]]-1))]
  xR = X[(Mphase[estar]:(Mphase[E[poskstar+1]]-1)),]
  y2 = array(c(yL,yR))
  x2 = rbind(xL,xR)
   
  ## Sample the edge vector to be maintained  for merge phase after removing the CP
  ## (newRight =1 if the edge vector to the Right of the CP is maintained, =0 otherwise)
  newRight = sample(0:1,1)
 
  if(nbVarMax>1){
    Sig2 = Sig2all[poskstar - 1 + newRight]
  } else {
    Sig2 = Sig2all
  }

  ## Update delta
  delta2 = sampleDelta2(poskstar-1+newRight, x2, q, Ball, Sall, Sig2, alphad2, betad2)

  ## Compute projection of the matrices required for the computation of the acceptation probability alpha
  #modif 17 avril 
  Px2 = computePx(length(y2), as.matrix(x2[,which(Sall[poskstar-1+newRight,] == 1)]), delta2) 
  PxL = computePx(length(yL), as.matrix(xL[,which(Sall[poskstar-1,] == 1)]), delta2)
  PxR = computePx(length(yR), as.matrix(xR[,which(Sall[poskstar,] == 1)]), delta2)

  ## Compute the acceptation probability alpha
  alpha = cp.computeAlpha(-1, sum(Sall[poskstar-1+newRight,])-1, s-1, Mphase[E[poskstar-1]], Mphase[estar], Mphase[E[poskstar+1]], yL, PxL, yR, PxR, y2, Px2, D, delta2, q, smax, v0, gamma0)
  
  ## Sample u to conclude either to  acceptation or to rejection
  u = runif(1,0,1)

  ## Boolean for the acceptation of the CP death move initially set to 0 (=1 if birth accepted, 0 otherwise)
  accept = 0
  
  if(!is.nan(alpha) & u <= alpha){
    ## Acceptation of the death of the selected CP
    ## Move acceptation boolean =1
    accept=1

    ## Position of the line of the phase to be removed in the matrices Sall and Ball
    away = (1-newRight) * poskstar + newRight * (poskstar-1)
   
    ## Remove the CP in E and the phase in the matrices Sall and Ball
    if(nbVarMax>1){
      Sig2all = Sig2all[(1:(s+1))[-c(away)]]
    }
    E = E[(1:(s+2))[-c(poskstar)]]
    Sall = matrix(Sall[(1:(s+1))[-c(away)],], s, q+1)
    Ball = matrix(Ball[(1:(s+1))[-c(away)],], s, q+1)
  }

  ##  Return all variables (+ variables 'move' and  accept describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases) and whether it was accepted  

  return( list( E=E, Sall=Sall, Ball=Ball, Sig2all=Sig2all, accept=accept, move=2, alpha=alpha, estar=estar))
}

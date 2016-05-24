cp.birth <-
function(E, Sall, Ball, Sig2all, X, Y, D, GLOBvar, HYPERvar){
  # INPUT: E, Sall, Ball ,Sig2all, X, Y, D, GLOBvar, HYPPERvar

  # current number of changepoints
  s = length(E) - 2
  
  ### assignement of global variables used here ###
  # ajoute par Sophie 18/09/09
  m=GLOBvar$m
  q = GLOBvar$q
  qmax = GLOBvar$qmax
  Mphase = GLOBvar$Mphase
  segMinLength = GLOBvar$segMinLength
  nbVarMax = GLOBvar$nbVarMax
  smax = GLOBvar$smax
  dyn = GLOBvar$dyn
  ### end assignement ###

  ### assignement of hyperparameters variables used here ###
  alphalbd = HYPERvar$alphalbd
  betalbd = HYPERvar$betalbd
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  ### end assignement ###
  

  ## search for possible CP, not in E and not close to E if segMinLength (length of phase) is > than 1
  toremove = E
   
  if(segMinLength>1) for(i in 1:(segMinLength-1)) toremove = c(toremove, E-i, E+i)
  # possible CPs are those not in 'toremove'
  possibleCP = setdiff((1+dyn):E[length(E)], toremove)
 
  
  ## Sample the new CP "estar"
  estar = sample(c(possibleCP, possibleCP),1)
  
  ## Position of the phase containing the new CP
  poskl = sum(E < estar)
  
  ## Current edges vector S in the phase containing the new CP
  Sold = Sall[poskl,]
  
  ## Current number of edges k in the phase containing the new CP
  k = sum(Sold) - 1

  ## Sample lambda
  lambda = rgamma(1, shape=alphalbd, rate=betalbd)

  ## Sample a new edges vector newS
  newS = array(1, q+1)
  newS[1:q] = 1:q %in% sample(1:q, sampleK(0, qmax, lambda, 1), replace=FALSE)

  ## Sample Right or Left (the position for the new edges vector: to the right or to the left of the new CP estar)
  if(runif(1,0,1) < 1/2){
    ##  New edges vector to the left of the new CP
    ## Boolean (= 1 if  the new edges vector is to the right of the new CP, 0 otherwise) 
    newRight = 0
    sL = newS
    sR = Sold
  } else {
    ## New edges vector to the right of the new CP
    ## Boolean (= 1 if  the new edges vector is to the right of the new CP, 0 otherwise)
    newRight = 1
    sR = newS
    sL = Sold
  }
  
  ## Compute the matrices required for the computation of the acceptation probability alpha
   
  yL = Y[(Mphase[E[poskl]]:(Mphase[estar]-1))]
  xL = X[(Mphase[E[poskl]]:(Mphase[estar]-1)),]
  yR = Y[(Mphase[estar]:(Mphase[E[poskl+1]]-1))]
  xR = X[(Mphase[estar]:(Mphase[E[poskl+1]]-1)),]
  y2 = array(c(yL, yR))
  x2 = rbind(xL, xR)
  
  ## Updating parameters
  # Variance Sig2
  Sig2 = Sig2all[poskl]

  # hyperparm delta2
  delta2 = sampleDelta2(poskl, x2, q, Ball, Sall, Sig2, alphad2, betad2)
  
  ## Compute projection of the matrices required for the computation of the acceptation probability alpha
  PxL = computePx(length(yL), as.matrix(xL[,which(sL == 1)]), delta2)
  PxR = computePx(length(yR), as.matrix(xR[,which(sR == 1)]), delta2)
  Px2 = computePx(length(y2), as.matrix(x2[,which(Sold == 1)]), delta2)
  
  ## Compute the acceptation probability alpha
  # modified by Sophie 18/09/09
  alpha = cp.computeAlpha(1, sum(newS)-1, s, Mphase[E[poskl]], Mphase[estar], Mphase[E[poskl+1]], yL, PxL, yR, PxR, y2, Px2, D, delta2, q, smax, v0, gamma0)

  
  ## Sample u to conclude either to  acceptation or to rejection
  u = runif(1,0,1)
  
  ## Boolean for the acceptation of the CP birth move initially set to 0 (=1 if birth accepted, 0 otherwise)
  accept = 0
  
  if(!is.nan(alpha) & u <= alpha){
    ## Acceptation of the birth of the new CP
    ## Move acceptation boolean =1
    accept=1
    
    ## Compute new Sig2
	newSig2 = array(0, s+2)
	newSig2[(1:(s+2))[-c(poskl,poskl+1)]] = Sig2all[(1:(s+1))[-c(poskl)]]
    
    ## Compute new regression parameters newB
    newB = matrix(0, s+2, q+1)
    newB[(1:(s+2))[-c(poskl,poskl+1)],] = Ball[(1:(s+1))[-c(poskl)],]
    
    ## Update newB 
    if(newRight == 0){
      ## Update the phase to the left of the new CP (in newB)
	  newSig2[poskl] = sampleSig2(yL,PxL,v0,gamma0)
	  newSig2[poskl+1] = Sig2all[poskl]
	  Sig2all = newSig2
	  Sig2 = newSig2[poskl]
      newB[poskl+1,] = Ball[poskl,]
      newB[poskl, which(newS == 1)] = sampleBxy(xL[,which(newS==1)], yL, Sig2, delta2)
    } else {
      ## Update the phase to the right of the new CP (in newB)
	  newSig2[poskl+1] = sampleSig2(yR, PxR, v0, gamma0)
	  newSig2[poskl] = Sig2all[poskl]
	  Sig2all = newSig2
	  Sig2 = newSig2[poskl]
      newB[poskl,] = Ball[poskl,]
      newB[poskl+1, which(newS == 1)] = sampleBxy(xR[, which(newS == 1)], yR, Sig2, delta2)
    }
    
    ## Update current model and parameters
    Ball = newB
    Sall = (abs(Ball)>0)*1
    E = sort(c(E,estar))
  }

  ##  Return all variables (+ variables 'accept' and 'move' describing the move type (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases) and whether it was accepted  
  return(list(E=E, Sall=Sall, Ball=Ball, Sig2all=Sig2all,  move=1, accept=accept, alpha=alpha, estar=estar))
}

bdu <-
function(u, rho3, x, y, S, Sig2, delta2, q, v0, gamma0){

  ### INPUT:u,rho,s=S[i,],sig2=Sig2[i],delta2.
  ###	x: the data in state i in columns 
  ###	ni: total nb of repeated measurements	
  ### OUTPUT: newS,newSig2,newB.
  ### depends on: 
  ### q the number of predictors
  ### constant v0, gamma0.

  ## Variable move describing the move type  (1= Edge birth, 2= Edge death, 3= Update coefficient, default=3)
  move = 3

  ## Boolean indicating whether the move is accepted or not (=1 if accepted, 0 otherwise, default=0)
  accept = 0

  ## New edges vector, to be returned at the end of the function
  newS = S


  
  ## Current number of edges
  l = sum(S) - 1 # = L[i]

  ## To update Sig2: matPx= Pxl, Pxlp1, or Pxlm1 depending the computed move (birth, death or update).
  ## Compute the projection matrix with the current edge ("Pxl")
  Pxl = computePx(length(y), x[,which(S == 1)], delta2)
 
  if(u < rho3[1]){
    ## Variable move describing the move type  (1= Edge birth, 2= Edge death, 3= Update coefficient)
    move = 1

    ## Sample the additional edge
    sstar = sample(c(which(S==0), which(S==0)), 1) # needed when there is only one position  S==0

    ## Proposed edges vector (with an additional edge)
    stmp = S
    stmp[sstar] = 1

    ## Compute the projection matrix with an additional edge ("Pxl plus 1")
    Pxlp1 = computePx(length(y), x[,which(stmp == 1)], delta2)

    ## Compute birth ratio
    rbirth = ((gamma0 + t(y) %*% Pxlp1 %*% y)/(gamma0 + t(y) %*% Pxl %*% y))^(-(length(y) + v0)/2)/sqrt(1 + delta2)

    ## Sample u 
    u = runif(1,0,1)
                                     
    if(u <= min(1,rbirth)){
      accept = 1
      newS = stmp
    }
    
    ## at this stage : 
    ## newS=S if birth rejected
    ##    =stmp if birth accepted
	
   } else {
    if(u < rho3[2]){
      ## Variable describing the move type  (1 for Edge birth, 2 for Edge death, 3 for Update coefficient)
      move=2
      
      ## Sample the added predictor
      sstar = sample(c(which(S[1:q]==1), which(S[1:q]==1)),1) # needed when there is only one position  S[1:q]==1
      ## Proposed edges vector (after taking away one edge)
      stmp = S
      stmp[sstar] = 0
      
      ## Compute the projection matrix after removimg one edge ("Pxl minus 1")
      Pxlm1 = computePx(length(y), x[,which(stmp==1)], delta2)

      ## Compute death ratio
      rdeath=((gamma0 + t(y) %*% Pxl %*% y)/(gamma0 + t(y) %*% Pxlm1 %*% y))^((length(y) + v0)/2)*(sqrt(1 + delta2))

      ## Sample u 
      u<-runif(1,0,1)
      
      if(u <= min(1,rdeath)){
        ## Boolean for the acceptation of the CP death move (=1 if birth accepted, 0 otherwise)
        accept = 1
        newS = stmp
       }
    }else{
      #updating coefficient only
      accept=1
    }
    
  }

 
  ## Updating coefficients (in all cases)
  newB = array(0, q+1)
  if(sum(newS) > 0){
    newB[which(newS == 1)] = sampleBxy(x[, which(newS==1)], y, Sig2, delta2)
  }
 
  ##  Return all variables
  return(list( newS=newS, newB=newB, u=u, move=move, accept=accept))
}

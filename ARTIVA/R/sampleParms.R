sampleParms <-
function(X, GLOBvar, HYPERvar, s=NULL, CPinit=NULL){
  ### assignement of global variables used here ###
  smax = GLOBvar$smax
  q = GLOBvar$q
  qmax = GLOBvar$qmax
  n = GLOBvar$n
  Mphase = GLOBvar$Mphase
  nbVarMax = GLOBvar$nbVarMax
  dyn = GLOBvar$dyn
  segMinLength = GLOBvar$segMinLength
  
  ### assignement of hyperparameters variables used here ###
  alphaD = HYPERvar$alphaD
  betaD = HYPERvar$betaD
  alphalbd = HYPERvar$alphalbd
  betalbd = HYPERvar$betalbd
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  ### end assignement ###

  ## Sample the number of breakpoint positions  
  if(!(is.null(CPinit))){
	E=CPinit
	s=length(E)-2
  }else{
	if( is.null(s) ){
		## If s=NULL, sample s.
		## sample D for the number of CP :
		D = rgamma(1, shape=alphaD, rate = betaD)  # scale s= 1/rate => f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)
    
		## Sample s
		s <- sampleK(0,smax,D,1)
	} else {
		## Update D 
		D = rgamma(1, shape=alphaD+s, rate = betaD+1)
	}

	## CP (phase > 2)
	E = c(1+dyn, n+1)
	cpt = s
	while(cpt > 0){
		# search for possible CP, not in E and not close to E if segMinLength (length of phase) is > than 1
		toremove = E
		if(segMinLength>1) for(i in 1:(segMinLength-1)) toremove = c(toremove, E-i, E+i)
		# possibles CPs are those not in 'toremove'
		possibleCP = setdiff((1+dyn):E[length(E)], toremove)
		# sample one CP in possibleCP (the vector is double for sake of function sample when size is = to 1)
                
  		if(length(possibleCP)==0){
                  cpt = 0
                  s = length(E)-2
                }else{
                  cp = sample( c(possibleCP, possibleCP), 1)              
                  E=sort(c(E, cp))
                  cpt = cpt-1
                }
	}
  }


  ###  sample model for each hidden state
  ## sample model structure
  S = matrix(0, s+1, q+1)
  for (i in 1:(s+1)){
    ## sample lambda
    lambda = rgamma(1, shape=alphalbd, rate = betalbd)  # scale s= 1/rate => f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)
    
    ## sample the nb of predictors
    kPred = sampleK(0, qmax, lambda, 1)

    if(kPred>0){
      S[i, sample(1:q, kPred, replace=FALSE)] = array(1, kPred) # structure du model (1 si pred in the model)
    }
  }

  ## we assume that there is a constante in each model
  S[, q+1] = array(1,s+1)

  ## sample sigma : IG(v0/2,gamma0/2)
  Sig2 = rinvgamma(n=min(nbVarMax,s+1), shape=v0/2, scale=gamma0/2)

  ## Coefficients 
  B = matrix(0,s+1,q+1)

  
  ## sample coef 
  for (i in 1:(s+1)){
    ## sample delta2
    delta2 = rinvgamma(1, shape=alphad2, scale=betad2)
     
    B[i,] = sampleBinit(S[i,], Sig2[i], delta2, X[(Mphase[E[i]]):(Mphase[E[i+1]]-1),], q)
    
  }
 
  return(list(E=E, S=S, B=B, Sig2=Sig2, s=s))
}

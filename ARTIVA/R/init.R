init <-
function(X, sinit, GLOBvar, HYPERvar, OUTvar, CPinit=NULL){
  ### assignement of global variables used here ###
  # storage matrix size = niter (maximal number of iterations before convergence) + PSRFactor_nbIterEval (nb iterations after convergence)
  niter = GLOBvar$niter
  if(OUTvar$PSRFactor)  niter = niter+OUTvar$PSRFactor_nbIterEval
    
  smax = GLOBvar$smax
  p = GLOBvar$p
  q = GLOBvar$q
  ### end assignement ###
  
  # stock along iteration
  # for stocking E (vector of CPs)
  Estock = matrix(0, niter, smax+p+1)
  # for stocking S (vector of network models -models coded in power2- for each possible phase (smax+1 phase possible))
  Sstock = matrix(-1, niter, smax+1)
  # for stocking B (vector of predictor coefficients + constant (q+1) for each phase (smax+1))
  Bstock = matrix(0, niter, (q+1)*(smax+p))
  # for stocking sigma (vector of sigma noise at each phase (smax+1))
  Sig2stock = matrix(0, niter, smax+1) 
  
  ## counting CP moves (CP Birth, CP death, CP move, Updating phases)
  CPMovesCount = array(0,4)
  # GL: Le 08-07-11: Names are added...
  names(CPMovesCount) = c("CPbirth", "CPdeath", "CPmove", "ModelUpdate")
  CPMovesAcceptation = array(0,4)
  names(CPMovesAcceptation) = c("CPbirth", "CPdeath", "CPmove", "ModelUpdate")


  ## counting "Updating phases" moves (Edge Birth, Edge death, Udating regression coefficient)
  EdgesMovesCount = array(0,3)
  names(EdgesMovesCount) = c("EdgeBirth", "EdgeDeath", "CoeffUpdate")
  EdgesMovesAcceptation = array(0,3)
  names(EdgesMovesAcceptation) = c("EdgeBirth", "EdgeDeath", "CoeffUpdate")
  counters = list(CPMovesCount=CPMovesCount, CPMovesAcceptation=CPMovesAcceptation, EdgesMovesCount=EdgesMovesCount, EdgesMovesAcceptation=EdgesMovesAcceptation)


  ## Initialisation
  init = sampleParms(X, GLOBvar, HYPERvar, s=sinit, CPinit) 

 
  ## Nb of changepoints
  s = init$s

  ## Vector of changepoints location
  E = init$E

  ## Boolean matrix (nb of phases x (q+1)) describing set of predictors of each phase
  Sall = init$S

  ## Real matrix (nb of phases x (q+1)) describing regression coefficients of each phase
  Ball = init$B

  ## Variance
  Sig2all = init$Sig2

  initState = list(s=s, E=E, Sall=Sall, Ball=Ball, Sig2all=Sig2all)
  
  ## Stock first iteration
  Estock[1,1:(s+2)] = init$E
  
  ##updated by Sophie 07/07/2011
  if(q==1){
    Sstock[1,1:(s+1)] = init$S[,1:q] 
  }else{
    Sstock[1,1:(s+1)] = init$S[,1:q] %*% 2^(0:(q-1))
  }
  
  Bstock[1,1:((s+1)*(q+1))] = t(init$B)	
  Sig2stock[1,1:length(init$Sig2)] = init$Sig2
  listStock = list(Estock=Estock, Sstock=Sstock, Bstock=Bstock, Sig2stock=Sig2stock)

  return(list(counters=counters, initState=initState, listStock=listStock))
}

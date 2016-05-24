#' Sample initial parameters for the MCMC simulation.
#' 
#' This function samples the initial hyperparameters and parameters that are
#' needed for the MCMC simulation.
#' 
#' 
#' @param X Input data.
#' @param GLOBvar Global variables of the MCMC simulation.
#' @param HYPERvar Hyperparameter variables.
#' @param s_init Initial number of changepoints.
#' @param options MCMC options, as given by e.g. \code{\link{defaultOptions}}.
#' @return Returns a list with elements: \item{E}{The initial changepoint
#' vector.} \item{S}{The intial networks structure.} \item{B}{The initial
#' regression parameters.} \item{Sig2}{The inital sigma squared variances.}
#' \item{betas}{The intial hyperparameters for the exponential information
#' sharing prior.} \item{hyper_params}{The initial hyperparameters for the
#' binomial information sharing prior.}
#' @author Sophie Lebre
#' 
#' Frank Dondelinger
#' @seealso \code{\link{init}}
#' @references For more information about the parameters and hyperparameters,
#' see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export sampleParms
sampleParms <-
function(X, GLOBvar, HYPERvar, s_init=NULL, options){
  ### Assignment of global variables used here ###
  smax = GLOBvar$smax
  q = GLOBvar$q
  qmax = GLOBvar$qmax
  n = GLOBvar$n
  Mphase = GLOBvar$Mphase
  nbVarMax = GLOBvar$nbVarMax
  dyn = GLOBvar$dyn
  minPhase = GLOBvar$minPhase
  beta_fixed = GLOBvar$beta_fixed
  ### End assignment ###
  
  ### Assignment of hyperparameter variables used here ###
  alphaD = HYPERvar$alphaD
  betaD = HYPERvar$betaD
  alphalbd = HYPERvar$alphalbd
  betalbd = HYPERvar$betalbd
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  ### End assignment ###

  E = list()
  s = matrix(0, q, 1)
 
  CPinit = options$cp.init 

  for(target in 1:q) {
    ## Sample the number of breakpoint positions  
    if(!(is.null(CPinit))){
      # Same CPs for all target variables
      if(!is.list(CPinit)) {
        E[[target]] = CPinit
      # Different CPs for different target variables
      } else {
        E[[target]] = CPinit[[target]]
      }
      
      s[target] = length(E[[target]]) - 2
    
    } else {
      
      if( is.null(s_init) ) {
        ## If s = NULL, sample s.
        ## Sample D for the number of CP :
        # scale s= 1/rate => f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)
        D = rgamma(1, shape=alphaD, rate=betaD)
        
        ## Sample s
        s[target] <- sampleK(0,smax,D,1)
        
      } else {
    	  ## Update D 
    	  D = rgamma(1, shape=alphaD+s, rate=betaD+1)
      }

      ## CP (phase > 2)
      E[[target]] = c(1+dyn, n+1)
      cpt = s[target]
    
      while(cpt > 0){
    	
        # Search for possible CP, not in E and not close to E 
        # if minPhase (length of phase) is > than 1
    	  toremove = E[[target]]
    	  
    	  if(minPhase>1) for(i in 1:(minPhase-1)) toremove = c(toremove, E[[target]]-i, E[[target]]+i)
    	  
        # Possibles CPs are those not in 'toremove'
    	  possibleCP = setdiff((1+dyn):E[[target]][length(E[[target]])], toremove)
    	  
        # Sample one CP in possibleCP (the vector is double for sake of 
        # function sample when size is = to 1)
    	  cp = sample( c(possibleCP, possibleCP), 1)
       
    	  E[[target]] = sort(c(E[[target]], cp))
    	  cpt = cpt-1
      }
    }
  }  
  
  S = list()
  B = list()
  Sig2 = list()
  
  # Intialise hyperparameters for information sharing
  if(is.null(options$hyper.init)) {
    betas = runif(q, 0, 1);
    hyper_params = matrix(1, 1, 4)
  } else {
    betas = options$hyper.init
    hyper_params = options$hyper.init
  }
  
  # Sample model structures from prior
  for(target in 1:length(X)) {
   
    S[[target]] = matrix(0, s[target]+1, q+1)
    
    for (i in 1:(s[target]+1)){
      ## Sample lambda
      lambda = rgamma(1, shape=alphalbd, rate = betalbd)  
      # scale s= 1/rate => f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)
    
      ## Sample the nb of predictors
      kPred = sampleK(0, qmax, lambda, 1)

      if(kPred>0){
        S[[target]][i, sample(1:q, kPred, replace=FALSE)] = array(1, kPred) 
        # structure of the model (1 if pred is in the model)
      }
      
      # Eliminate self-loops if not allowed
      if(!GLOBvar$self.loops) {
        S[[target]][i, target] = 0
      }
      
      # Set fixed edges
      if(!is.null(GLOBvar$fixed.edges)) {
        fixed.edges = GLOBvar$fixed.edges
        fixed.indices = c(fixed.edges[,target] != -1)
        S[[target]][i,c(fixed.indices, FALSE)] = 
          fixed.edges[fixed.indices,target]
      }
    }
    
    ## We assume that there is a constant in each model
    S[[target]][, q+1] = array(1,s[target]+1)
    
    # Sample hyperparameters
  
    ## Sample sigma : IG(v0/2,gamma0/2)
    Sig2[[target]] = rinvgamma(n=min(nbVarMax,s[target]+1), shape=v0/2, 
                               scale=gamma0/2)
    
    ## Regression Coefficients 
    B[[target]] = matrix(0,s[target]+1,q+1)

    ## Sample coefficients
    for (i in 1:(s[target]+1)){
      if(nbVarMax == 1){ iSig = 1 } else { iSig = i }
      
      ## sample delta2
      delta2 = rinvgamma(1, shape=alphad2, scale=betad2)

      B[[target]][i,] = sampleBinit(S[[target]][i,], 
              Sig2[[target]][iSig], delta2, 
              X[[target]][(Mphase[E[[target]][i]]):(Mphase[E[[target]][i+1]]-1),], 
              q)
    }
  }
  
  return(list(E=E, S=S, B=B, Sig2=Sig2, s=s, betas=betas, 
              hyper_params=hyper_params))
}


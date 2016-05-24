#' Initialise the MCMC simulation.
#' 
#' This function intialises the parameters and variables needed for the MCMC
#' simulation.
#' 
#' 
#' @param X Input response data.
#' @param Y Input target data.
#' @param sinit Initial changepoints.
#' @param GLOBvar Global variables used during the MCMC simulation.
#' @param HYPERvar Hyperparameter variables.
#' @param options MCMC simulation options as obtained e.g. by
#' \code{\link{defaultOptions}}.
#' @return List with elements: \item{counters}{Matrices for counting the number
#' of moves made and accepted.} \item{initState}{Initial state of the variables
#' of the MCMC simulation.} \item{listStock}{Variables for recording the
#' network, changepoint and hyperparameter samples.}
#' @author Sophie Lebre
#' 
#' Frank Dondelinger
#' @seealso \code{\link{sampleParms}}
#' @export init
init <-
function(X, Y, sinit, GLOBvar, HYPERvar, options){
  ### Assignment of global variables used here ###
  niter = GLOBvar$niter
  smax = GLOBvar$smax
  p = GLOBvar$p
  q = GLOBvar$q
  birth_proposals = GLOBvar$birth_proposals 
  method = GLOBvar$method
  prior.params = GLOBvar$prior.params
  ## end assignment ###
  
  # Record of CPs
  niter = min(niter, 10000)
  
  hyperstock = matrix(0, niter, q)
  Bstock = list(); Estock = list()
  
  ## Counting moves 
  # CP Birth, CP death, CP move, Edge Move, Hyperparameter Level 1 Move, 
  # Hyperparameter Level 2 Move
  cptMove = array(0,6)
  acceptMove = array(0,6)
  
  counters = list(cptMove=cptMove, 
                  acceptMove=acceptMove)
  
  for(target in 1:length(X)) {
    # For stocking B (vector of predictor coefficients + 
    # constant (q+1) for each phase (smax+1))
    Bstock[[target]] = matrix(0, niter-1, (q+1)*(smax+p))
    # For stocking E (vector of changepoints)
    Estock[[target]] = matrix(0, niter-1, smax+p+1)
  }

  ## Initialisation
  init = sampleParms(X, GLOBvar, HYPERvar, s_init=sinit, options) 

  ## Nb of changepoints
  s = init$s

  ## Vector of changepoints location
  E = init$E

  # Boolean matrix (nb of phases x (q+1)) describing 
  # the set of predictors of each phase
  Sall = init$S
 
  # Real matrix (nb of phases x (q+1)) describing 
  # the regression coefficients of each phase
  Ball = init$B
  
  k = 0
  
  # Hyperparameters (influence of hypernetwork)
  if(is.null(prior.params)) { 
      
    # Exponential Model
    if(method == 'exp_soft' || method == 'exp_hard') {
      prior.params = init$betas;
    # Binomial Model
    } else if(method == "bino_hard" || method == "bino_soft") {
      prior.params = init$hyper_params
    # No information sharing model
    } else {
      prior.params = -1
    }
  }
  
  # Proposal width (only used for exponential IS models)
  hyper.proposals = 0.01
  
  ## Variance
  Sig2all = init$Sig2

  initState = list(s=s, E=E, Sall=Sall, Ball=Ball, Sig2all=Sig2all, 
                   prior.params=prior.params, k=k, 
                   hyper.proposals=hyper.proposals)  

  ## Stock first iteration
  for(target in 1:length(X)) {
    Estock[[target]][1,1:(s[target]+2)] = init$E[[target]]
    Bstock[[target]][1,1:((s[target]+1)*(q+1))] = c(t(init$B[[target]]))
  }
  
  listStock = list(Estock=Estock, Bstock=Bstock, hyperstock=hyperstock)
  
  return(list(counters=counters, initState=initState, listStock=listStock))
}


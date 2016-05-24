#' Main function of the MCMC simulation.
#' 
#' This function executes the main loop of the MCMC simulation, making the
#' different moves and recording samples.
#' 
#' 
#' @param X Input response data.
#' @param Y Input target data.
#' @param initiation Initialisation of the MCMC simulation, as obtained by
#' function \code{\link{init}}.
#' @param GLOBvar Global variables of the MCMC simulation.
#' @param HYPERvar Hyperparameter variables.
#' @return Returns a list with the following elements: \item{counters}{List
#' containing the different move counters for the number of times moves have
#' been proposed and accepted.} \item{listStock}{List containing the recorded
#' samples for the networks, changepoints and hyperparameters}
#' @author Sophie Lebre
#' 
#' Frank Dondelinger
#' @seealso \code{\link{runDBN}}
#' @references For more information about the MCMC simulations, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export main
main <-
function(X, Y, initiation, GLOBvar, HYPERvar){

  ### Assignment of global variables used here ###
  niter = GLOBvar$niter
  smax = GLOBvar$smax
  q = GLOBvar$q
  birth_proposals = GLOBvar$birth_proposals
  sequential_model = GLOBvar$sequential_model
  method = GLOBvar$method
  ### End assignment ###

  ### Assignment of hyperparameters variables used here ###
  cD = HYPERvar$cD
  alphaD = HYPERvar$alphaD
  betaD = HYPERvar$betaD
  ### End assignment ###
  
  ### Assignment of initiation variables used here ###
  # Initial state
  E = initiation$initState$E
  Sall = initiation$initState$Sall
  Ball = initiation$initState$Ball
  Sig2all = initiation$initState$Sig2all
  s = initiation$initState$s
  
  # Counters
  cptMove = initiation$counters$cptMove
  acceptMove = initiation$counters$acceptMove

  # Storage matrices
  Estock = initiation$listStock$Estock
  Sstock = initiation$listStock$Sstock
  Bstock = initiation$listStock$Bstock
  hyperstock = initiation$listStock$hyperstock;
  Sig2stock = initiation$listStock$Sig2stock
  
  # Hyperparameters for information sharing prior
  HYPERvar$prior.params = initiation$initState$prior.params
  HYPERvar$hyper.proposals = initiation$initState$hyper.proposals
  ### End assignment ###
  
  # How often to monitor the acceptance rate for hyperparameters proposal
  # tuning
  monitorRate = 500

  samples = 2:niter
  
  deltastock = matrix(0, length(samples), 1);
  chi2stock = hyperstock;
  
  i = 1;
  
  if(dim(Estock[[1]])[1] < niter) {
    samples = sort(sample(2:niter, dim(Estock[[1]])[1]))
    samples[dim(Estock[[1]])[1]] = niter
  }
  
  
  # Do niter interations
  for (r in 2:niter){
    
	  ## Print percentage of accomplished iteration 
    if ((r %% (round(niter/20))) == 0) {
	    cat(round(5 * r/(round(niter/20))), "% ")
        
      if(GLOBvar$method != 'poisson')
        print(round(HYPERvar$prior.params, digits=2))
      
      if(i > 5000 && GLOBvar$psrf.check) {
        psrf_r = c()
        
        for(target in 1:q) {
          psrf_r = c(psrf_r, psrf_check(Bstock[[target]], q, smax, i)) 
        }
        
        print(psrf_r)
        
        print(psrf_check_hyper(hyperstock, i))
      }
      
    }
      
    target = sample(1:q, 1)
    
    D = rgamma(1, shape=s[target]+alphaD, rate = 1+betaD)

    rho = computeRho4(s[target], 0, smax, cD, D)

    # Sample u to choose one of the 4 moves : 
    # CP birth, CP death, CP shift, Update phases.
    u1 = runif(1, 0, 1)

    # Choose 1 out of the 4 moves (depending on the value of u)
    if (u1 < rho[1]){
      # CP birth move: return the new model if the move 
      # is accepted, the previous model otherwise.
      out = cp.birth(E, Sall, Ball, Sig2all[[target]], 
                     X[[target]], Y[[target]], D, GLOBvar, 
                     HYPERvar, target)
    } else if(u1 < rho[2]){
      # CP death move: return the new model if the move is accepted, 
      # the previous model otherwise.
      out = cp.death(E, Sall, Ball, Sig2all[[target]], X[[target]], 
                     Y[[target]], D, GLOBvar, HYPERvar, target)
    } else if(u1 < rho[3]){
      # CP shift move: return the new model if the move is accepted, 
      # the previous model otherwise.
      out = cp.shift(E, Sall[[target]], Ball[[target]], Sig2all[[target]], 
                     X[[target]], Y[[target]], GLOBvar, HYPERvar, target)
    } else {
      # Update phases: return the new model if the move is accepted, 
      # the previous model otherwise. (includes hyperparameter moves)
        
      out = phase.update(E, Sall, Ball, Sig2all, X, Y, GLOBvar, 
          HYPERvar, target)
    }

    ## Apply changes to the current model
    E[[target]] = out$E
    
    if(is.list(out$Ball)) {
      Ball = out$Ball
    } else {
      Ball[[target]] = out$Ball
    }
    
    for(node in 1:q) {
      Sall[[node]] = (abs(Ball[[node]])>0)*1
      s[node] = length(E[[node]])-2
    }

    HYPERvar$prior.params = out$prior.params;
    HYPERvar$k            = out$k
    
    if(r %% monitorRate == 0) {
      HYPERvar$hyper.proposals = 
        proposalTuning(acceptMove/cptMove, HYPERvar$hyper.proposals)
    }
    
    if(is.list(out$Sig2all)) {
      Sig2all = out$Sig2all  
    } else {
      Sig2all[[target]] = out$Sig2all
    }
    
    ## Update moves counts
    cptMove[out$move] = cptMove[out$move]+1
    acceptMove[out$move] = acceptMove[out$move]+out$accept

    ## Stock model and parameters
    if(r == samples[i]) {      
      
      for(node in 1:q) {
        Estock[[node]][i,1:(s[node]+2)] = E[[node]]
        Bstock[[node]][i,1:((s[node]+1)*(q+1))] = c(t(Ball[[node]]))
      }
      
      hyperstock[i,1:length(HYPERvar$prior.params)] = out$prior.params

      i = i + 1;
    }
    
    
    
  } # end iteration

  counters = list(cptMove=cptMove, 
                  acceptMove=acceptMove)
  listStock = list(Estock=Estock, Bstock=Bstock, 
                   hyperstock=hyperstock, deltastock=deltastock, chi2stock=chi2stock)

  return(list(counters=counters, listStock=listStock))
}


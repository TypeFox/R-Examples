MNS <-
function(dat, lambda_pop, lambda_random, parallel=FALSE, cores=NULL, max_iter=100, tol=1e-5){
  # estimate a population graphical model using MIXED NEIGHBOURHOOD SELECTION
  #
  # The approach here builds on the neighbourhood selection method introduced by N Meinshausen & P Buhlmann 2006.
  # However we account for the heterogenious nature of data by using penalized linear mixed models as opposed to 
  # Lasso models which assume data is IID
  # This allows us to understand both the average "population" relationships between nodes but also where the
  # variability is concentrated across subjects
  #
  #
  # INPUT:
  #      - dat: a list where each entry is the data for a given subject  
  #      - lambda_pop & lambda_random: penalty parameters for fixed and random effects respectively
  #        Note that if we setp lambda_pop=0 then we are only penalizing the random effects!
  #      - parallel: boolean indicating whether to estimate model in parallel. Default is false
  #        Similarly, cores indicates the number of cores to use (if null the default value from doParallel is selected)
  #      - max_iter, tol: max number of iterations & convergence tolerance
  #
  # OUTPUT:
  #      - PresPop: population neighbourhood selection network
  #      - PresRE: estimate of edge specific random effects
  #      - PresBLUP: 
  #
  #
  #
  
  ## ADD SAFETY CHECK HERE
  
  # define some variables:
  p = ncol(dat[[1]]) # number of variables
  n = nrow(dat[[1]]) # number of observations for each subject (assumed fixed across subjects)
  N = length(dat) # number of subjects
  subject = rep(1:N, each=n) # subject encoding - not actually used anywhere but may be used in future!
  mIter = 0 # mean iteration over nodes
  NLL = vector("list", p)
  
  # put data in correct format:
  o = PrepareData(dat) 
  # prepare output:
  PresPop = matrix(0, ncol=p, nrow=p) # population network
  PresRE = matrix(0, ncol=p, nrow=p) # estimate of the noise of each edge
  PresBLUP = array(0, c(p,p,N)) # BLUPs for each subject
  PresNoise = matrix(0, ncol=p, nrow=p)
  
  if (parallel){
    # we proceed to apply in parallel using plyr package
    if (is.null(cores)) cores=2
    
    # make parallel instance:
    workers=makeCluster(2)
    registerDoParallel(workers,cores=cores)
    
    #manually export things I will need
    clusterExport(cl=workers, varlist =list("PenalizedLMM_opt", "ginv", "BuildRandomDesign", "glmnet", "calculateNLL", "dmvnorm", "subject", "lambda_pop", "lambda_random"),envir=environment())
    
    #results = llply(o, .fun = function(xentry){
    #  y = xentry[[1]]
    #  X = xentry[[2]]
    #  subNeighbourhood = PenalizedLMM_opt(y = y, X = X, subject = subject, lambda_pop = lambda_pop, lambda_random = lambda_random)
    #  return(subNeighbourhood)
    #},  .parallel=TRUE)
    
    results = parLapply(cl = workers, X = o, fun = function(xentry){
       y = xentry[[1]]
       X = xentry[[2]]
       subNeighbourhood = PenalizedLMM_opt(y = y, X = X, subject = subject, lambda_pop = lambda_pop, lambda_random = lambda_random)
       return(subNeighbourhood)
      })
    
    stopCluster(workers)
    # now tidy up results and put into correct format:

    # store results:
    for (node in 1:p){
      PresPop[node, -node]= results[[node]]$beta
      PresRE[node, -node] = results[[node]]$sigmaRE
      PresNoise[node, -node] = results[[node]]$sigmaNoise
     
      for (i in 1:N){
        PresBLUP[,,i][node, -node] = results[[node]]$BLUPs[i,]
      }
      
      mIter = mIter + results[[node]]$it
      NLL[[node]] = results[[node]]$NLL
    }
    
  } else {
    # now we proceed to estimate the neighbourhood of each node SEQUENTIALLY
    for (node in 1:p){
      y = o[[node]][[1]]
      X = o[[node]][[2]]
      #Z = o[[node]][[3]]
      #subNeighbourhood = PenalizedLMM_opt(y = y, X = X, Z = Z, subject = subject, lambda_pop = lambda_pop, lambda_random = lambda_random)
      subNeighbourhood = PenalizedLMM_opt(y = y, X = X, subject = subject, lambda_pop = lambda_pop, lambda_random = lambda_random)
      
      # store results:
      PresPop[node, -node]= subNeighbourhood$beta
      PresRE[node, -node] = subNeighbourhood$sigmaRE
      PresNoise[node, -node] = subNeighbourhood$sigmaNoise
      
      for (i in 1:N){
        PresBLUP[,,i][node, -node] = subNeighbourhood$BLUPs[i,]
      }
      
      mIter = mIter + subNeighbourhood$it
      NLL[[node]] = subNeighbourhood$NLL
    }
  }
  
  MNSobj = list(PresPop=PresPop, PresRE=PresRE, PresBLUP=PresBLUP, PresNoise=PresNoise, it=mIter/p, NLL=NLL)
  class(MNSobj) = "MNS"
  
  return(MNSobj)
}

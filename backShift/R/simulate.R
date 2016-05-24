#' Simulate data of a causal cyclic model under shift interventions.
#' 
#' @description Simulate data of a causal cyclic model under shift interventions.
#' 
#' @param n Number of observations.
#' @param p Number of variables.
#' @param A Connectivity matrix A. The entry \eqn{A_{ij}} contains the edge 
#'        from node i to node j.
#' @param G Number of environments, has to be larger than two for \code{backShift}.
#' @param intervMultiplier Regulates the strength of the interventions.
#' @param noiseMult Regulates the noise variance.
#' @param nonGauss Set to \code{TRUE} to generate non-Gaussian noise.
#' @param hiddenVars Set to \code{TRUE} to include hidden variables.
#' @param knownInterventions Set to \code{TRUE} if location of interventions
#'  should be known.
#' @param fracVarInt If \code{knownInterventions} is \code{TRUE}, fraction of 
#'  variables that are intervened on in each environment.
#' @param simulateObs If \code{TRUE}, also generate observational data.
#' @param seed Random seed.
#' @return A list with the following elements: 
#' \itemize{
#'   \item \code{X} (nxp)-dimensional data matrix
#'   \item \code{environment} Indicator of the experiment or the intervention type an 
#'   observation belongs to. A numeric vector of length n. 
#'   \item \code{interventionVar} (Gxp)-dimensional matrix with intervention variances.
#'   \item \code{interventions} Location of interventions if \code{knownInterventions}
#'    was set to \code{TRUE}.
#'   \item \code{configs} A list with the following elements: 
#'   \itemize{
#'   \item \code{trueA} True connectivity matrix used to generate the data.
#'   \item \code{G} Number of environments.
#'   \item \code{indexObservationalData} Index of observational data
#'   \item \code{intervMultiplier} Multiplier steering the intervention strength
#'   \item \code{noiseMult} Multiplier steering the noise level
#'   \item \code{fracVarInt}  If \code{knownInterventions} was set to \code{TRUE}, 
#'    fraction of variables that were intervened on in each environment.
#'   \item \code{hiddenVars} If \code{TRUE}, hidden variables exist.
#'   \item \code{knownInterventions} If \code{TRUE}, location of interventions is known.
#'   \item \code{simulateObs} If \code{TRUE}, environment \code{1} contains 
#'    observational data.
#'   }
#' }
#' @references Dominik Rothenhaeusler, Christina Heinze, Jonas Peters, Nicolai Meinshausen (2015):
#' backShift: Learning causal cyclic graphs from unknown shift interventions.
#' arXiv preprint: \url{http://arxiv.org/abs/1506.02494}
simulateInterventions <- function(n, p, A, G, intervMultiplier, noiseMult, 
                                  nonGauss, hiddenVars, knownInterventions, 
                                  fracVarInt, simulateObs,  seed =1){
  ####### set seed
  set.seed(seed)
  
  ####### initialize
  X <- Perturb <- matrix(0,nrow=n,ncol=p)
  envVar <- matrix(0, nrow = G, ncol = p)  
  if(simulateObs) idxObs <- 1 else idxObs <- NULL
  
  ###### simulate
  if(!knownInterventions){
    ## intervention targets are unknown
    environment <- rep(1:G, each=ceiling(n/G))[1:n]
    
    # additive soft interventions
    
    ## simulate noise perturbations in each environment
    for (i in unique(environment)){
      ind <- which(environment==i)
      if(!(simulateObs && i == idxObs)){
        # if observational data should be simulated: no interventions in environment 1
        multiplier <- rexp(p)*intervMultiplier
        envVar[i,] <- multiplier
        Perturb[ind,] <- sweep(matrix(rnorm(length(ind)*p),ncol=p), 2, multiplier,FUN="*")
      }
    }
  
    interventions <- NULL
  }else{
    
    ## choose explicity intervention targets
    unique.interventions <- list()
    for (g in 1:G){
      unique.interventions[[g]] <- sample(1:p,round(p*fracVarInt))
      envVar[g, unique.interventions[[g]]] <- 
        rexp(length(unique.interventions[[g]]))*intervMultiplier
    }
    
    interventions <- list()
    for (i in 1:n) interventions[[i]] <- sample(unique.interventions,1)[[1]]
    environment <- match(interventions, unique(interventions))
    
    if(simulateObs){
      obsDataPoints <- which(environment == idxObs) 
      for(i in seq_along(obsDataPoints))
        interventions[[ obsDataPoints[i] ]] <- numeric(0)
    }
    
    # additive soft interventions
    ## change level of noise for each intervention
    for (i in 1:n){
      if(is.null(idxObs) || !(simulateObs & environment[i] == idxObs)){
        # if observational data should be simulated: no interventions in environment 1
        Perturb[i, interventions[[i]]] <- 
          rnorm(length(interventions[[i]]))*envVar[environment[i],interventions[[i]]]
      }
    }
    
  }

  if(hiddenVars){
    ## Input of hidden variables into each variable
    gamma <- rnorm(p)*noiseMult
    W <- if(nonGauss){
      expo <- matrix(rexp(n),nrow=n)
      rs <- matrix(1, nrow = n)
      rs[sample(n, size = 0.5*n)] <- -1
      rs * expo
    }else{
      rnorm(n)
    }
    Input <- matrix(outer(W, gamma, FUN="*"), nrow = n)
  }else{
    ## Independent noise at each variable 
    Input <- if(nonGauss){
        expo <- matrix(rexp(n*p),nrow=n)
        rs <- matrix(1, nrow = n, ncol = p)
        rs[sample(n*p, size = 0.5*n*p)] <- -1
        rs * expo
      }else{
        matrix(rnorm(n*p),nrow=n)
      }
  }  
  
  inv <- solve(diag(p) - A)
  X <- (Input + Perturb)%*%inv
  
  config.options <- list(trueA = A, G = G, indexObservationalData = idxObs,
                         intervMultiplier = intervMultiplier, 
                         noiseMult = noiseMult, fracVarInt = fracVarInt, 
                         hiddenVars = hiddenVars, 
                         knownInterventions = knownInterventions, 
                         simulateObs = simulateObs)
  
  list(X = X, environment = environment, 
       interventionVar = envVar^2, 
       interventions = interventions, 
       configs = config.options)
}

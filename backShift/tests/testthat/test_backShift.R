library(backShift)
context("All supported arguments")

# Example --------

seed <- 1
set.seed(seed)

## Parameters for simulation ------

# number of observations
n <- 100
# number of variables
p <- 10
# number of environments
G <- 5
# if the location of the interventions is known, on how many vars. should
# be intervention in each environment (as a fraction of p)
fracVarInt <- 0.5
# multiplier for interventions (m_I in manuscript)
intMult <- 1.5
# multiplier for interventions (m_e in manuscript)
noiseMult <- 1

data("exampleAdjacencyMatrix")
A <- exampleAdjacencyMatrix
p <- 10

## Options for method -------

# if stability selection should not be used, set EV = 0
# number of false selections for stability selection
EV <- 2
# selection threshold for stability selection
thres <- 0.75

boolean <- c(TRUE, FALSE)

# simulate non-Gaussian noise? 
for(nonGauss in boolean){
  # also simulate observational data?
  for(simulateObs in boolean){
    # should hidden vars be included?
    for(hidden in boolean){
      # should the location of the interventions be known?
      for(knownInterventions in boolean){
          ## Simulate data -------
          simulation.res <- simulateInterventions(n, p, A, G, intMult, noiseMult, 
                                                  nonGauss = nonGauss, 
                                                  hiddenVars = hidden,
                                                  knownInterventions = knownInterventions, 
                                                  fracVarInt, 
                                                  simulateObs = simulateObs, 
                                                  seed)
          
          
          # extract X, environment vector and index of observational data
          X <- simulation.res$X
          env <- simulation.res$environment
          baseInd <- simulation.res$configs$indexObservationalData
          
          # use covariance matrix instead of Gram matrix
          for(useCov in boolean){
            
          ## Run backShift -------
          test_that(paste("Checks backShift for.."), {
            expect_is(
              backshift.res <- backShift(X, env, covariance=useCov, nsim = 5,
                                         ev=EV, threshold=thres, 
                                         baseSettingEnv = baseInd, tolerance = 1e-3, 
                                         verbose = FALSE)$Ahat
              , "matrix")
            
          })
          
        }
        
        
      }
    }
  }
}
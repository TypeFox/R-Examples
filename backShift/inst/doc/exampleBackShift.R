## ---- message = FALSE----------------------------------------------------
# load package
require(backShift)
require(ggplot2)
require(fields)

# set seed
seed <- 1
set.seed(seed)

## ------------------------------------------------------------------------
# use covariance matrix instead of Gram matrix
useCov <- TRUE 

## ------------------------------------------------------------------------
# bound on expected number of false selections for stability selection
EV <- 2
# selection threshold for stability selection
thres <- 0.75 

## ------------------------------------------------------------------------
# threshold for point estimate 
thres.pe <- 0.25

## ---- results = "hide"---------------------------------------------------
# number of variables
p <- 10

# set whether to provide or to generate the adjacency matrix A 
providedA <- TRUE

if(providedA){
  data("exampleAdjacencyMatrix")
  A <- exampleAdjacencyMatrix
  p <- 10
}else{
  # parameters to generate A
  
  # should A be cyclic?
  cyclic <- TRUE
  # expected number of neighbors per node
  expNumNeigh <- 0.1*p 
  # range for coefficients
  minCoef <- 0.3
  maxCoef <- 0.8
  
  ## Generate A -------
  cat("Generating A...\n") 
  A.gen.result <- generateA(p, expNumNeigh, minCoef, maxCoef, cyclic)
  A <- A.gen.result$A
  cat("A has a cycle of size", A.gen.result$sizeCycle, "\n") 
}

## ---- results = "hide"---------------------------------------------------
# number of observations
n <- 10000
# number of environments
G <- 10 
# also simulate observational data?
simulateObs <- TRUE 
# should hidden vars be included?
hidden <- FALSE 
# should the location of the interventions be known?
knownInterventions <- FALSE 
# if the location of the interventions is known, how many vars. should
# be intervened on in each environment (as a fraction of p)
fracVarInt <- 0.5
# multiplier for interventions (m_I in manuscript)
intMult <- 1.5
# multiplier for interventions (m_e in manuscript)
noiseMult <- 1
# simulate non-Gaussian noise? 
nonGauss <- FALSE

## Simulate data -------
cat("Simulating the data...\n") 
simulation.res <- simulateInterventions(n, p, A, G, intMult, noiseMult, 
                                        nonGauss, hidden, knownInterventions, 
                                        fracVarInt, simulateObs, seed)

# extract X, environment vector and index of observational data
X <- simulation.res$X
env <- simulation.res$environment
baseInd <- simulation.res$configs$indexObservationalData

## ------------------------------------------------------------------------
## Run backShift -------
backshift.res <- backShift(X, env, covariance=useCov, ev=EV, threshold=thres, 
                        baseSettingEnv = baseInd, tolerance = 1e-6, 
                        verbose = FALSE)

## ---- fig.show='hold', fig.width=6, fig.height=6, fig.align='center', results = "hide"----
# extract estimates
Ahat <- backshift.res$Ahat
Ahat.structure <- backshift.res$AhatAdjacency

# compute performance metrics and plot result
cat("Plotting true graph... \n") 
plotGraphEdgeAttr(estimate = A, plotStabSelec = FALSE, labels = colnames(A), 
                  thres.point = 0, thres.stab = thres, main = "True graph")

## ---- fig.show='hold', fig.width=6, fig.height=6, fig.align='center', results="hide"----
cat("Plotting point estimate, thresholded at", thres.pe,"... \n") 
plotGraphEdgeAttr(estimate = Ahat, plotStabSelec = FALSE, labels = colnames(A), 
                  thres.point = thres.pe, thres.stab = thres, 
                  main = paste("Point estimate thresholded at", thres.pe))

## ---- fig.show='hold', fig.width=6, fig.height=6, fig.align='center', results = "hide"----
cat("Plotting stability selection result... \n") 
plotGraphEdgeAttr(estimate = Ahat.structure, plotStabSelec = TRUE, 
                  labels = colnames(A), thres.point = thres.pe, 
                  edgeWeights = Ahat, thres.stab = thres, 
                  main = "Stability selection result")

## ---- results = "asis"---------------------------------------------------
# metrics for point estimate, thresholded at thres.pe
metricsThresholdedA <- metricsThreshold(A, Ahat, thres = thres.pe)

# metrics for stability selection result
metricsStabSelection <- metricsThreshold(A, Ahat.structure, thres = 0)

## ---- results = "asis", echo = FALSE-------------------------------------
colnames(metricsThresholdedA)[3] <- colnames(metricsStabSelection)[3] <- "TPR/Recall"
knitr::kable(metricsThresholdedA, digits = 4, caption = paste("Metrics for point estimate, thresholded at", thres.pe))
knitr::kable(metricsStabSelection, digits = 4, caption = "Metrics for stability selection result")

## ---- message=FALSE, fig.width=10, fig.height=7.5------------------------
plotInterventionVars(backshift.res$varianceEnv, simulation.res$interventionVar)

## ---- fig.keep='all', fig.width=3, fig.height=3--------------------------
for(i in 1:G){
  plotDiagonalization(estConnectivity = backshift.res$Ahat, X = X, env = env, whichEnv = i)
}


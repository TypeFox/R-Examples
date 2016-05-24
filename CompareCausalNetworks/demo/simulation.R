## 2nd example --------
# Simulate data with connectivity matrix A with assumptions
# 1) No hidden variables
# 2) Precise location of interventions is known

if(!requireNamespace("backShift", quietly = TRUE))
 stop("The package 'backShift' is needed for the examples to 
 work. Please install it.", call. = FALSE)

require(backShift)

## simulate data
myseed <- 1

# sample size n
n <- 10000

# p=5 predictor variables and connectivity matrix A
p <- 5
labels <- c("1", "2", "3", "4", "5")
A <- diag(p)*0
A[1,2] <- 0.8
A[2,3] <- -0.8
A[3,4] <- 0.8
A[3,5] <- 0.8
A[4,5] <- 0.3

# can add/remove feedback by using/not using
# A[5,2] <- 0.8 

# divide data in 10 different environments
G <- 10

# simulate choose explicity intervention targets
simResult <- simulateInterventions(n, p, A, G, intervMultiplier = 3, 
             noiseMult = 1, nonGauss = TRUE, hiddenVars = FALSE, 
             knownInterventions = TRUE, fracVarInt = 0.2, 
             simulateObs = TRUE, seed = myseed)
X <- simResult$X
environment <- simResult$environment
interventions <- simResult$interventions

# number of unique environments
G <- length(unique(environment))

## apply all  methods given in vector 'methods'
methods <- c("ICP", "hiddenICP", "gies")

# select whether you want to run stability selection
stability <- FALSE

# arrange graphical output into a rectangular grid
sq <- ceiling(sqrt(length(methods)+1))
par(mfrow=c(ceiling((length(methods)+1)/sq),sq))

## plot and print true graph
cat("\n true graph is  ------  \n" )
print(A)
plotGraphEdgeAttr(A, plotStabSelec = FALSE, labels = labels, 
 thres.point = 0, main = "TRUE GRAPH")

## loop over all methods and compute and print/plot estimate
for (method in methods){
  cat("\n result for method", method,"  ------  \n" )
 
  if(!stability){
    # Option 1): use this estimator as a point estimate
    Ahat <- getParents(X, environment, interventions = interventions,
                       method=method, alpha=0.1, pointConf = TRUE)
  }else{
    # Option 2): use a stability selection based estimator
    # with expected number of false positives bounded by EV=2
    Ahat <- getParentsStable(X, environment,  interventions = interventions,
                             EV=2, method=method, alpha=0.1)
  }
 
  # print and plot estimate (point estimate thresholded if numerical estimates
  # are returned)
  print(Ahat)
  if(!stability)
    plotGraphEdgeAttr(Ahat, plotStabSelec = FALSE, labels = labels, 
     thres.point = 0.05,
     main=paste("POINT ESTIMATE FOR METHOD\n", toupper(method)))
  else
    plotGraphEdgeAttr(Ahat, plotStabSelec = TRUE, labels = labels,
     thres.point = 0, main = paste("STABILITY SELECTION ESTIMATE\n FOR METHOD", 
     toupper(method)))
 }
context('structure-moves')

dataset = simulateNetwork(l=30, cps=c(10,20))

# Prepare for structure move

num.iter = 500
  
sim.data = dataset$sim_data
options = defaultOptions()
options$cp.fixed = TRUE
options$cp.init = list()
options$cp.init[[1]] = c(2,7,17,dim(sim.data)[2])
options$cp.init[[2]] = c(2,13,23,dim(sim.data)[2])
  
for(i in 3:dim(sim.data)[1]) {
  options$cp.init[[i]] = c(2,dim(sim.data)[2])
}
  
information.sharing='poisson'
prior.params=NULL
  
data = sim.data
  
# Time series length
n = dim(data)[2]
  
# Number of variables
q = dim(data)[1]
  
##### Important remark ####
# if you know the changepoint position, and you want to run the procedure only for estimating the model within phases, 
# you can  set the variable 'CPinit' to the known CP position  set the variable 'cD=0' in the file 'hyperparms.R' 
# then CP move will never be considered (only the 4th move 'phase.update' will be considered).
# however there is still some probleme with the 'output' fonctions in this case, I can help updating this when I will be back to work. 
  
targetdata = data; preddata=NULL; multipleVar = TRUE; minPhase=2;
niter = num.iter; scaling = TRUE; method = information.sharing; 
self.loops = TRUE; k = 15; outputFile='.';
  
nbCPinit=min(floor(n/2),5)
  
m = options$m
dyn = options$dyn
  
# Position of each time point in the data (designed for the algorithm)
Mphase = seq(1,n*m+1,by=m)-dyn*m
  
if(multipleVar){
  nbVarMax=options$maxCP+1
}else{ 
  nbVarMax=1
}
  
# Read input data
targetData = targetdata
predData = targetData
  
# Standardise inputs to N(0,1)
if(scaling){
  targetData = t(scale(t(targetData)))
  predData = t(scale(t(predData)))
}
  
fixed.edges = matrix(-1, q, q)

# A few tests :
# targetData and predData must have the same number of columns
if(ncol(targetData) != ncol(predData)) stop("Target data and input data don't have the same number of columns.\n")
# The number of columns corresponds to n (timepoints) x m (repetitions)
if(ncol(targetData) != n*m) stop("Number of columns incompatible with n and m.\n")
  
# List of genes analyzed :
# Analyze all rows of targetData
posResponse = 1:nrow(targetData)
  
# Names of predictors
# Take rownames of predData
predNames = row.names(predData)
  
# Names of targets
# Take rownames of predData
targetNames = row.names(targetData)
  
# Position of the predictor variables in the data for each response 
# (matrix [nrow(targetData) x q]) 
# By default all the predictors of predData are taken for each gene
bestPosMat = matrix(1:q, nrow(targetData), q, byrow=TRUE)
  
### Create Global Variables used in all functions
GLOBvar = list(n=n, m=m, p=1, q=q, qmax=options$maxTF, smax=options$maxCP, 
                 dyn=options$dyn, 
                 minPhase=minPhase, nbVarMax=nbVarMax, Mphase=Mphase, bestPosMat=bestPosMat, 
                 niter=niter, target=NULL,predNames=predNames, targetNames=targetNames, 
                 lmax=options$lmax, method=method, 
                 prior.params=prior.params, self.loops=self.loops,
                 burnin=options$burnin, psrf.check=options$psrf.check,
                 pp.l1=options$pp.l1, pp.l2=options$pp.l2, hyper.fixed=options$hyper.fixed,
                 cp.fixed=options$cp.fixed, fixed.edges=fixed.edges)
  
### Create HyperParms Variables used in all functions
HYPERvar = HyperParms(options)
HYPERvar$k = k
  
### Create Output Variables used in output functions
OUTvar = list(outputFile=outputFile, by.node=options$save.by.node,
                save.file=options$save.file)
  
X = list()
Y = list()
  
# For each target variable, prepare input
for(target in posResponse) {
  GLOBvar$target=target
    
  ## Build response variables Y and predictor variables X
  input = buildXY(targetData, predData, GLOBvar)
    
  X[[target]] = input$X
  Y[[target]] = input$Y   
}
  
## Initialize system
initiation = init(X, Y, nbCPinit, GLOBvar, HYPERvar, options)
  
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
Eall = initiation$initState$E
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


## Test function for rejection changes
noChangeOnReject <- function() {
  change.flag = TRUE
  

  
  # Do niter interations
  for (r in 2:niter){
    
    target = sample(1:q, 1)
    
    D = rgamma(1, shape=smax+alphaD, rate = 1+betaD)
    
    rho = computeRho4(s[target], 0, smax, cD, D)
    
    # Sample u to choose one of the 4 moves : 
    # CP birth, CP death, CP shift, Update phases.
    u1 = runif(1, 0, 1)
    
    E = Eall[[target]]
    # Current number of changepoints
    s = length(E) - 2
    
    ### Assignment of global variables used here ###
    q = GLOBvar$q
    qmax = GLOBvar$qmax
    Mphase = GLOBvar$Mphase
    nbVarMax = GLOBvar$nbVarMax
    smax = GLOBvar$smax
    lmax = GLOBvar$lmax
    method = GLOBvar$method
    small_prop = GLOBvar$small_prop
    self.loops = GLOBvar$self.loops
    ### End assignment ###
    
    ### Assignment of hyperparameters variables used here ###
    c = HYPERvar$c
    alphalbd = HYPERvar$alphalbd
    betalbd = HYPERvar$betalbd
    alphad2 = HYPERvar$alphad2
    betad2 = HYPERvar$betad2
    v0 = HYPERvar$v0
    gamma0 = HYPERvar$gamma0
    k = HYPERvar$k
    ### End assignment ###
    
    prior.params = HYPERvar$prior.params
    
    Sig2 = Sig2all[[target]]
    
    ## Observations in the chosen phase
    y = Y[[target]]
    x = X[[target]]
    
    model = 0
    
    # Group information about network segments    
    network.info = CollectNetworkInfo(Sall, Eall, prior.params, -1, 
                                      target, q, self.loops, k)
    

    ## Compute a structure move 
    bduout =  make_structure_move(x, y, Sall[[target]], Ball[[target]], Sig2, q, 
                                    qmax, network.info, method, Mphase, E, 
                                  fixed.edges[,target], HYPERvar)
      
    change.flag = change.flag && 
      (bduout$accept || all(bduout$newS == Sall[[target]]))
    
  } # end iteration
  
 
  return(change.flag)
}

## Test function for updating network information
updateCorrectly <- function() {
  same.flag = TRUE
  
  target = sample(1:q, 1)  
    
  E = Eall[[target]]
  # Current number of changepoints
  s = length(E) - 2
    
  ### Assignment of global variables used here ###
  q = GLOBvar$q
  self.loops = GLOBvar$self.loops
    ### End assignment ###
    
  k = HYPERvar$k
    
  prior.params = HYPERvar$prior.params
    
  # Group information about network segments    
  network.info = CollectNetworkInfo(Sall, Eall, prior.params, -1, 
                                      target, q, self.loops, k)

  S.new = Sall[[target]]
  S.new[sample(1:dim(S.new)[1],1),sample(1:dim(S.new)[2], 1)] = (runif(1) > 0.5)*1
  Sall[[target]] = S.new[,, drop=FALSE]

  network.info = addProposalNetworkInfo(network.info, S.new[,1:q,drop=FALSE], E)
  
  network.info.new = CollectNetworkInfo(Sall, Eall, prior.params, -1, 
                                          target, q, self.loops, k)
  
  for(segment.i in 1:length(network.info$nets)) {
    same.flag = same.flag && 
      all(network.info.new$nets[[segment.i]] == network.info$new.nets[[segment.i]])
  }
  
  return(same.flag)
}

## Test function for mapping from local to global segments in network information
getMapping <- function() {
  same.flag = TRUE
  
  target = sample(1:q, 1)  
  
  E = Eall[[target]]
  # Current number of changepoints
  s = length(E) - 2
  
  ### Assignment of global variables used here ###
  q = GLOBvar$q
  self.loops = GLOBvar$self.loops
  ### End assignment ###
  
  k = HYPERvar$k
  
  prior.params = HYPERvar$prior.params
  
  # Group information about network segments    
  network.info = CollectNetworkInfo(Sall, Eall, prior.params, -1, 
                                    target, q, self.loops, k)
  
  return(network.info$global.mapping[2,])
}

# Test for mapping from local to global segments within network.info
test_that('mapping from local to global segments correct',
          expect_that(getMapping(), equals(c(1,1,2,2,3))))

# Test for updating network info within structure move
test_that('network info the same before and after ',
          expect_that(updateCorrectly(), is_true()))

# Test for potential bug with local and global segments
test_that('rejected moves make no change',
          expect_that(noChangeOnReject(), is_true()))

# Test output is not NULL (Generic Run Test) 
test_that('output not null',
          expect_that(is.null(EDISON.run(sim.data, num.iter=num.iter)), is_false()))

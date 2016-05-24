#' Setup and run the MCMC simulation.
#' 
#' This function initialises the variabes for the MCMC simulation, runs the
#' simulation and returns the output.
#' 
#' 
#' @param targetdata Target input data: A matrix of dimensions NumNodes by
#' NumTimePoints.
#' @param preddata Optional: Input response data, if different from the target
#' data.
#' @param q Number of nodes.
#' @param n Number of timepoints.
#' @param multipleVar \code{TRUE} when a specific variance is estimated for
#' each segment, \code{FALSE} otherwise.
#' @param minPhase Minimal segment length.
#' @param niter Number of MCMC iterations.
#' @param scaling If \code{TRUE}, scale the input data to mean 0 and standard
#' deviation 1, else leave it unchanged.
#' @param method Network structure prior to use: \code{'poisson'} for a sparse
#' Poisson prior (no information sharing), \code{'exp_hard'} or
#' \code{'exp_soft'} for the exponential information sharing prior with hard or
#' soft node coupling, \code{'bino_hard'} or \code{'bino_soft'} with hard or
#' soft node coupling.
#' @param prior.params Initial hyperparameters for the information sharing
#' prior.
#' @param self.loops If \code{TRUE}, allow self-loops in the network, if
#' \code{FALSE}, disallow self-loops.
#' @param k Initial value for the level-2 hyperparameter of the exponential
#' information sharing prior.
#' @param options MCMC options as obtained e.g. by the function
#' \code{\link{defaultOptions}}.
#' @param outputFile File where the output of the MCMC simulation should be
#' saved.
#' @param fixed.edges Matrix of size NumNodes by NumNodes, with 
#' \code{fixed.edges[i,j]==1|0} if the edge between nodes i and j is fixed, and 
#' -1 otherwise. Defaults to \code{NULL} (no edges fixed).
#' @return A list containing the results of the MCMC simulation: network
#' samples, changepoint samples and hyperparameter samples. For details, see
#' \code{\link{output}}.
#' @author Sophie Lebre
#' 
#' Frank Dondelinger
#' @seealso \code{\link{output}}
#' @references For more information about the MCMC simulations, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export runDBN
runDBN <-
function(targetdata, preddata=NULL, q, n,
  multipleVar=TRUE, minPhase=2, 
  niter=20000, scaling=TRUE, 
  method='poisson', prior.params=NULL,
  self.loops=TRUE, k = 15, options=NULL, outputFile='.',
  fixed.edges=NULL) {
  # runtvDBN: run whole program
  # Description:
  #
  # Arguments:
  # targetdata = target data (either the name of a file, or directly a matrix)
  # preddata = optional, file with predictor data when differing from the target data (either the name of a file, or directly a matrix), default=NULL.
  # q = number of predictor variables
  # n = number of timepoints in the series (do not mention de repetitions)
  # m = number of repeated measurements for each timepoint
  # p = number of response variables (default=1)
  # dyn = delay (in timepoints) between TF and target (default=0)
  # multipleVar = TRUE when a specific variance is estimated for each phase, FALSE otherwise (default: TRUE).
  # minPhase = minimal length of a phase (>1 if no repetition, default: 2)
  # maxCP = maximal number of CP (default= min(n-1,15))
  # maxTF = maximal number of TF (default=min(q,15))
  # nbCPinit = number of CP at initialisation (default=min(floor(n/2),5))
  # CPinit = optional, initial chengpoint vector, default=NULL 
  # alphaCP, betaCP = hyperparms for sampling the number k of CP : k ~ Gamma(alphaCP,betaCP)  (default: alphaCP=0.5, betaCP =1). You can use function choosePriors to set alphaCP and  betaCP according to the desired dimension penalisation.
  # alphaTF, betaTF= hyperparms for sampling the number l of TF : l ~ Gamma(alphaTF,betaTF) (default: alphaTF=0.5, betaTF =1).  You can use function choosePriors to set alphaTF and  betaTF according to the desired dimension penalisation.
  # pkCP, pkTF = prior distribution for the nulmber of CP or TF (default=NULL, necessary when BFOut = TRUE and the hyperparameters alpha, beta are note among the one provided by the package, see choosePriors for the list of available hyperparameters).
  # posResponse = row position (in targetdata) of targets to be analyzed (default: analize all genes)
  # bestPosMatFile = file containing row position of predictors for each gene (see default below)
  # niter = number of iterations (default: 20000)
  # ndeb, nfin = start and end of the iterations considered in the analyzes (default: 1 and niter)
  # simpleOut = TRUE for simple output  (the chosen model is the most represented model)
  # BFOut = TRUE for Bayes Factor (BF) analysis (selected model has the number of changepoints with the highest BF and, for eacg pahse the number of TF with the highest BF) 
  # Picture = do you want output pictures ?
  # Format = output picture format (1=jpg, 0=eps)
  # WriteStock = do you want to store all results (all values of all parameters over the iterations) ?
  # secondModel = do you want the results of the second best model ?
  # predNamesFile = file containing the names of the predictor (matrix format) (by default rownames of preddata will be used)
  # outputFile = name of output file (default: tvDBNoutput)
  
  nbCPinit=min(floor(n/2),5)
  
  if(is.null(options)) {
    options = defaultOptions()
  }
  
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
	predData=targetData

  # Standardise inputs to N(0,1)
  if(scaling){
  	targetData = t(scale(t(targetData)))
  	predData = t(scale(t(predData)))
  }
  
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
  print("Initialisation successful.")

  ## Run niter iterations
  print("Starting tvDBN iterations...")
  runiteration = main(X, Y, initiation, GLOBvar, HYPERvar)
  print("")
  print("End of iterations")
    
    
  ## Save output data
  result = output(runiteration$counters, runiteration$listStock, GLOBvar, HYPERvar, OUTvar)

  return(result)
}


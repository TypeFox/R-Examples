############################################################
# ReferenceClass class that holds problem on master process
############################################################

krigeProblem <- setRefClass("krigeProblem",
                            
                            fields = list(
                              localProblemName = "character",
                              
                              n = "numeric",
                              m = "numeric",
                              h_n = "numeric",
                              h_m = "numeric",
                              
                              meanFunction = "function",
                              predMeanFunction = "function",
                              covFunction = "function",
                              crossCovFunction = "function",
                              predCovFunction = "function",
                              
                              data = "ANY",  # by setting to ANY, this allows the value to be NULL allows so that can be left NULL
                              params = "ANY", # by setting to ANY, this allows the value to be NULL allows so that can be left NULL
                              
                              meanCurrent = "logical",
                              predMeanCurrent = "logical",
                              postMeanCurrent = "logical",
                              covCurrent = "logical",
                              crossCovCurrent = "logical",
                              predCovCurrent = "logical",
                              postCovCurrent = "logical",
                              cholCurrent = "logical",
                              predCholCurrent = "logical",
                              postCholCurrent = "logical",
                              
                              logDens = "numeric",
                              
                              optimOutput = "list",
                              inputs = "list"),
                            
                            methods = list(
                              initialize = function(localProblemName = NULL, numProcesses = NULL, h_n = NULL, h_m = NULL,  n = length(data), m = NULL, meanFunction = function(){}, predMeanFunction = function(){}, covFunction = function(){}, crossCovFunction = function(){}, predCovFunction = function(){}, inputs = NULL, params = NULL, data = NULL, packages = NULL, parallelRNGpkg = "rlecuyer", seed = 0, ...){
                                'Description:
Constructor that initializes the krigeProblem object. The krigeProblem class includes methods for standard kriging calculations and metadata necessary for carrying out the methods in a distributed fashion.


Arguments:

localProblemName: character. A string to be used to refer to the distributedKrigeProblem objects on the slave processes.

numProcesses: numeric. The number of processes amongst which to distribute the objects and calculations. Should be equal to D*(D+1)/2 for some integer D.

h_n: numeric. Block replication factor for the observation locations. If NULL, a good default value will be used.

h_m: numeric. Block replication factor for the prediction locations. If NULL, a good default value will be used.

n: numeric. Number of observation locations. Should match the length of "data".

m: numeric. Number of prediction locations.

meanFunction: function. A function that calculates the mean at the observation locations. See help on krigeProblem for more details on the requirements for this function.

predMeanFunction: function. A function that calculates the mean at the prediction locations. See help on krigeProblem for more details on the requirements for this function.

covFunction: function. A function that calculates the covariance matrix for the observation locations. See help on krigeProblem for more details on the requirements for this function.

crossCovFunction: function. A function that calculates the cross-covariance matrix between the observation and prediction locations. See help on krigeProblem for more details on the requirements for this function.
 
predCovFunction: function. A function that calculates the covariance matrix for the prediction. See help on krigeProblem for more details on the requirements for this function.

inputs: list. Contains components used in the functions above, used in conjunction with the indices contained on each process.

params: numeric. Vector (possibly named) of parameters.

data: numeric. Vector of values at the observation locations.

packages: character. Character vector with the names of packages that need to be loaded on the slave processes for the krigeProblem mean and covariance functions.

parallelRNGpkg: character. Character string indicating the parallel random number generation package to use. Should be either "rlecuyer" (the recommended default), "rsprng", or NULL. The latter is not guaranteed to prevent overlap between the random number streams on the separate processes, used in simulating realizations. 
 
seed: numeric. Seed used for parallel random number generation on the processes. Using the same seed should give one the same realizations, allowing for reproducibility.

...: other arguments to be passed to "setClass".


Details:
This constructor initializes the krigeProblem object on the master process and associated distributedKrigeProblem objects on the slave processes. It calculates the mean vector and covariance matrix for the observation locations based on the parameters. 


Value:
krigeProblem. The initialized object.
'
                                
                                if(is.null(localProblemName))
                                  stop("krigeProblem$new: You must specify a name, 'localProblemName', for the problem on the workers.")
                                
                                localProblemName <<- localProblemName
                                
                                bigGP.init(numProcesses, parallelRNGpkg = parallelRNGpkg, seed = seed)
                                
                                if(is.null(n)){
                                  stop("krigeProblem$new: You must specify a value for 'n'.")
                                } else n <<- n
                                if(is.null(m)){
                                  m <<- 0
                                  rm(m)
                                } else m <<- m
                                
                                if(is.null(h_n)){
                                  h_n <<- calcH(n)
                                  rm(h_n) # remove local variable so object's field is used below
                                } else h_n <<- h_n  
                                
                                if(is.null(h_m)){
                                  if(m != 0) h_m <<- calcH(m) else h_m <<- 1
                                  rm(h_m)
                                } else h_m <<- h_m
                                
                                meanFunction <<- meanFunction
                                predMeanFunction <<- predMeanFunction
                                covFunction <<- covFunction
                                crossCovFunction <<- crossCovFunction
                                predCovFunction <<- predCovFunction
                                inputs <<- inputs
                                initializeSlaveProblems(packages)
                                
                                data <<- data
                                if(!is.null(data)){
                                  if(length(data) != n) stop("krigeProblem$new: Length of 'data' should be equal to 'n'.")
                                  distributeVector(data, objPos = localProblemName, n = n, h = h_n)
                                }
                                
                                if(!is.null(params)){
                                  setParams(params)
                                  remoteConstructMean(pred = FALSE)
                                  remoteConstructCov(pred = FALSE, cross = FALSE)
                                
                                  cat("Initial mean and covariance construction complete.\n")
                                }

                                
                                callSuper(...)
                              },
                              
                              calcH = function(n){
                                'Description:
Calculates a good choice of the block replication factor given the number of values.


Arguments:

n: numeric. The number of values being subdivided amongst the processes and blocks.


Value:
numeric. The value of H, the block replication factor.
'
                                return(max(n %/% (.bigGP$D * 1000), 1))
                              },
                              
                              show = function(verbose = TRUE){
                                'Description:
Show method for krigeProblem objects, reporting dimensions of the problem, processes, and blocking, as well as status of the various distributed objects.
'
                                cat("bigGP container for a distributed kriging problem,\nstored on worker processes as '", localProblemName, "',\nwith a partition factor (D) of ", .bigGP$D, " and ", .bigGP$P, " worker processes (P):\n\n", sep = "")
                                cat(n, " observed locations (n) and ", ifelse(is.null(m), 0, m), " prediction locations (m).\n", sep = "")
                                cat("Block replication factors are ", h_n, ", ", ifelse(is.null(h_m), 0, h_m), " for the observation (h_n) and prediction (h_m) dimensions.\n",sep = "")
                                if(verbose){
                                  cat("Mean vector is", ifelse(meanCurrent, " ", " NOT "), "current.\n", sep = "")
                                  cat("Prediction (prior) mean vector is", ifelse(predMeanCurrent, " ", " NOT "), "current.\n", sep = "")
                                  cat("Posterior mean vector is", ifelse(postMeanCurrent, " ", " NOT "), "current.\n", sep = "")
                                  cat("Covariance matrix is", ifelse(covCurrent, " ", " NOT "), "current.\n", sep = "")
                                  cat("Crosscovariance matrix is", ifelse(crossCovCurrent, " ", " NOT "), "current.\n", sep = "")
                                  cat("Posterior covariance matrix is", ifelse(postCovCurrent, " ", " NOT "), "current.\n", sep = "")
                                  cat("Cholesky is", ifelse(cholCurrent, " ", " NOT "), "current.\n", sep = "")
                                  cat("Prediction (prior) cholesky is", ifelse(predCholCurrent, " ", " NOT "), "current.\n", sep = "")
                                  cat("Posterior cholesky is", ifelse(postCholCurrent, " ", " NOT "), "current.\n", sep = "")
                                }
                                invisible(.self)
                              },
                              
                              
                              initializeSlaveProblems = function(packages){
                                'Description:
Sets up the slave processes to carry out the krigeProblem distributed calculations.


Arguments:

packages: character vector. Names of packages that need to be loaded on the slave processes in order for the krigeProblem functions (meanFunction, predMeanFunction, covFunction, crossCovFunction, predCovFunction) to run.

Value:
none
'
                                for(pkg in packages) {
                                  status <- mpi.remote.exec(require, pkg, ret = TRUE)
                                  if(sum(unlist(status)) != length(status))
                                    stop(paste("krigeProblem$initialize: error in loading package ", pkg, " on slaves.", sep = ""))
                                }
                                tmp <- list(h_n = h_n, h_m = h_m, n = n, m = m, meanFunction = meanFunction, predMeanFunction = predMeanFunction, covFunction = covFunction, crossCovFunction = crossCovFunction, predCovFunction = predCovFunction, inputs = inputs, params = params)
                                mpi.bcast.Robj2slave(tmp)
                                mpi.bcast.cmd(tmp2 <- distributedKrigeProblem$new(tmp$h_n, tmp$h_m, tmp$n, tmp$m, tmp$meanFunction, tmp$predMeanFunction, tmp$covFunction, tmp$crossCovFunction, tmp$predCovFunction, tmp$inputs, tmp$params))
                                status <- mpi.remote.exec(localAssign, localProblemName, 'tmp2', ret = TRUE)
                                if("try-error" %in% sapply(status, class))
                                  stop("localAssign: error on slaves:\n", status)
                                
                                remoteRm(tmp, tmp2)
                                invisible(NULL)
                              },
                              
                              setParams = function(params, verbose = TRUE){
                                'Description:
Set parameter values in krigeProblem.


Arguments:

params: numeric. Vector of parameter values in form (possibly a named vector) for use in the krigeProblem mean and covariance functions.

verbose: logical. Should status information be printed?


Value:
none
'
                                if(!is.numeric(params))
                                  stop("krigeProblem$setParams: 'params' must be numeric.")
                                params <<- params
                                meanCurrent <<- FALSE
                                covCurrent <<- FALSE
                                cholCurrent <<- FALSE
                                predMeanCurrent <<- FALSE
                                postMeanCurrent <<- FALSE
                                crossCovCurrent <<- FALSE
                                predCovCurrent <<- FALSE
                                postCovCurrent <<- FALSE
                                predCholCurrent <<- FALSE
                                postCholCurrent <<- FALSE
                                
                                push(params, objPos = localProblemName)
                                if(verbose)
                                  cat("Set parameter values.\n")
                                invisible(NULL)
                              },
                              
                              remoteConstructMean = function(obs = TRUE, pred = !obs, verbose = FALSE){
                                'Description:
Calculate mean vector(s) in a distributed fashion.


Arguments:

obs: logical. Should mean be calculated for observed locations?

pred: logical. Should mean be calculated for prediction locations?

verbose: logical. Should status information be printed?


Value:
none
'
                                
                                if(is.null(params))
                                  warning("krigeProblem$remoteConstructMean: There are no parameter values.")
                                status <- mpi.remote.exec(localKrigeProblemConstructMean, localProblemName, obs = obs, pred = pred, ret = TRUE)
                                if("try-error" %in% sapply(status, class))
                                  stop("localKrigeProblemConstructMean: error on slaves:\n", status)

                                if(verbose)
                                  cat("Constructed requested mean vector(s).\n")
                                invisible(NULL)
                              },
                              
                              
                              remoteConstructCov = function(obs = TRUE, pred = FALSE, cross = FALSE, verbose = FALSE){
                                'Description:
Calculate covariance matrix/matrices in a distributed fashion.


Arguments:

obs: logical. Should covariance matrix be calculated for observed locations?

pred: logical. Should covariance matrix be calculated for prediction locations?

cross: logical. Should cross-covariance matrix between observed and prediction locations be calculated?

verbose: logical. Should status information be printed?


Value:
none
'
                                if(is.null(params))
                                  stop("krigeProblem$remoteConstructCov: There are no parameter values.")
                                status <- mpi.remote.exec(localKrigeProblemConstructCov, localProblemName, obs = obs, pred = pred, cross = cross, ret = TRUE)
                                if("try-error" %in% sapply(status, class))
                                  stop("localKrigeProblemConstructCov: error on slaves:\n", status)

                                if(verbose)
                                  cat("Constructed requested covariance matrix.\n")
                                invisible(NULL)
                              },
                              
                              
                              calcLogDeterminant = function( ) {
                                'Description:
Calculate log determinant of covariance matrix for observed locations.


Arguments:
none


Value:
numeric. The value of the log determinant.
'
                                if(!cholCurrent){
                                  remoteCalcChol(matName = 'C', cholName = 'L', matPos = localProblemName, cholPos = localProblemName, n = n, h = h)
                                  cholCurrent <<- TRUE
                                }
                                diagVals <- collectDiagonal("L", localProblemName, n = n, h = h_n)
                                return(2 * sum(log(diagVals)))
                              },
                              
                              
                              calcLogDens = function(newParams = NULL, newData = NULL, negative = FALSE, verbose = TRUE){
                                'Description:
Calculate log density of a set of values according to the Gaussian process.


Arguments:

newParams: numeric. Vector of parameter values. If NULL, parameter values stored as "params" in the krigeProblem object are used.

newData: numeric. Vector of data values. If NULL, data values stored as "data" in the krigeProblem object are used.

negative: logical. Should the negative log density be calculated?

verbose: logical. Should status information be printed?


Value:
numeric. The value of the log density.
'
                                if(!is.null(newParams))
                                  setParams(newParams, verbose)
                                if(!is.null(newData))
                                  distributeVector(newData, objName = 'data', objPos = localProblemName, n = n, h = h_n)
                                
                                if(!meanCurrent){
                                  remoteConstructMean(obs = TRUE, pred = FALSE, verbose = verbose)
                                  meanCurrent <<- TRUE
                                }
                                if(!covCurrent){
                                  remoteConstructCov(obs = TRUE, pred = FALSE, cross = FALSE, verbose = verbose)
                                  covCurrent <<- TRUE
                                }
                                if(!cholCurrent){
                                  remoteCalcChol(matName = 'C', cholName = 'L', matPos = localProblemName, cholPos = localProblemName, n = n, h = h_n)
                                  cholCurrent <<- TRUE
                                }
                                remoteCalc("data", "mean", `-`, 'tmp', input1Pos = localProblemName, input2Pos = localProblemName, outputPos = ".GlobalEnv")
                                remoteForwardsolve(cholName = 'L', inputName = 'tmp', outputName = 'result', cholPos = localProblemName, n1 = n, h1 = h_n)
                                ssq <- sum( collectVector('result', n = n, h = h_n)^2 )
                                logDens <<- -0.5 * (calcLogDeterminant() + ssq)
                                if(verbose)
                                  if(negative){
                                    cat("Negative log density value: ", -logDens, "; Parameters: ", params, ".\n")
                                  } else{
                                    cat("Log density value:", logDens, "; Parameters:", params, ".\n", sep = " ")
                                  }
                                if(negative) return(-logDens) else return(logDens)
                              },
                              
                              
                              optimizeLogDens = function(newParams = NULL, newData = NULL, method = "Nelder-Mead", verbose = FALSE, gr = NULL, lower = -Inf, upper = Inf, control = list(), hessian = FALSE, ...) {
                                'Description:
Optimizes the log density with respect to the parameter values.


Arguments:

newParams: numeric. Vector of parameter values. If NULL, parameter values stored as "params" in the krigeProblem object are used.

newData: numeric. Vector of data values. If NULL, data values stored as "data" in the krigeProblem object are used.

method: character string. The optimization method to be used. See help for optim().

gr: function. Gradient function for the log density. See help for optim().

lower, upper: numeric. Bounds on the variables for the "L-BFGS-B" method, or bounds in which to _search_ for method "Brent". See help for optim().

control: list. A list of control parameters for the optimizationGradient function for the log density. See help for optim().

hessian: logical. Should a numerically differentiated Hessian matrix be returned? See help for optim().

...: Further arguments to be passed to the function being optimized and to "gr". See help for optim().


Value:
list. Components are those retured by optim().
'
                                if(!is.null(newParams))
                                  setParams(newParams, verbose)
                                if(!is.null(newData))
                                  distributeVector(newData, objName = 'data', objPos = localProblemName, n = n, h = h_n)
                                
                                if(is.null(params))
                                  stop("krigeProblem$optimizeLogDens: There are no parameter values.")
                                optimOutput <<- optim(params, .self$calcLogDens, gr = gr, negative = TRUE, verbose = verbose, ..., method = method, lower = lower, upper = upper, control = control, hessian = hessian) # not clear why I need .self here
                                if(optimOutput$convergence == 0){
                                  calcLogDens(optimOutput$par) # this sets params and recalcs density as optim may not have used final parameter values in last call to calcLogDens
                                  cat("Optimization result (see optimOutput field for result of call to optim()):\nMLE: ", params, "\nObjective: ", optimOutput$value, "\n")
                                } else{
                                  cat("Convergence failure: see optimOutput field for details\n")
                                }
                                invisible(optimOutput)
                              },
                              
                              
                              predict = function(ret = FALSE, se.fit = FALSE, verbose = FALSE){
                                'Description:
Makes predictions (posterior mean values) at prediction locations, and optionally estimates prediction standard errors (i.e., posterior standard deviations).


Arguments:

ret: logical. Should the predictions be returned? They may not be needed if prediction is an interim calculation, with values distributed amongst processes to be used for later calculations.

se.fit: logical. Should prediction standard errors be calculated?

verbose: logical. Should status information be printed?


Value:
numeric. The vector of predictions, or NULL if "ret" is FALSE.
'
                                
                                if(!postMeanCurrent){
                                  if(!meanCurrent){
                                    remoteConstructMean(obs = TRUE, pred = FALSE, verbose = verbose)
                                    meanCurrent <<- TRUE
                                  }
                                  if(!predMeanCurrent){
                                    remoteConstructMean(obs = FALSE, pred = TRUE, verbose = verbose)
                                    predMeanCurrent <<- TRUE
                                  }
                                  if(!covCurrent){
                                    remoteConstructCov(obs = TRUE, pred = FALSE, cross = FALSE, verbose = verbose)
                                    covCurrent <<- TRUE
                                  }
                                  if(!cholCurrent){
                                    remoteCalcChol(matName = 'C', cholName = 'L', matPos = localProblemName, cholPos = localProblemName, n = n, h = h_n)
                                    cholCurrent <<- TRUE
                                  }
                                  remoteCalc("data", "mean", `-`, 'tmp1', input1Pos = localProblemName, input2Pos = localProblemName)
                                  remoteForwardsolve(cholName = 'L', inputName = 'tmp1', outputName = 'tmp2', cholPos = localProblemName, n1 = n, h1 = h_n)
                                  remoteBacksolve(cholName = 'L', inputName = 'tmp2', outputName = 'tmp3', cholPos = localProblemName, n1 = n, h1 = h_n)
                                  if(!crossCovCurrent){
                                    remoteConstructCov(obs = FALSE, pred = FALSE, cross = TRUE, verbose = verbose)
                                    crossCovCurrent <<- TRUE
                                  }
                                  remoteCrossProdMatVec(matName = 'crossC', inputName = 'tmp3', outputName = 'tmp4', matPos = localProblemName, n1 = n, n2 = m, h1 = h_n, h2 = h_m)
                                  remoteCalc("predMean", 'tmp4', `+`, "postMean", input1Pos = localProblemName, outputPos = localProblemName)
                                  postMeanCurrent <<- TRUE
                                  if(verbose)
                                    cat("Calculated predicted values (i.e., posterior mean at prediction locations).\n")
                                } else 
                                  if(verbose)
                                    cat("Predicted values are current; returning values from child processes if requested.\n")
                                if(se.fit){
                                  if(ret){
                                    if(!covCurrent){
                                      remoteConstructCov(obs = TRUE, pred = FALSE, cross = FALSE, verbose = verbose)
                                      covCurrent <<- TRUE
                                    }
                                    if(!cholCurrent){
                                      remoteCalcChol(matName = 'C', cholName = 'L', matPos = localProblemName, cholPos = localProblemName, n = n, h = h_n)
                                      cholCurrent <<- TRUE
                                    }
                                    if(!crossCovCurrent){
                                      remoteConstructCov(obs = FALSE, pred = FALSE, cross = TRUE, verbose = verbose)
                                      crossCovCurrent <<- TRUE
                                    }
                                    remoteForwardsolve(cholName = 'L', inputName = 'crossC', outputName = 'tmp1', cholPos = localProblemName, inputPos = localProblemName, n1 = n, n2 = m, h1 = h_n, h2 = h_m)
                                    remoteCrossProdMatSelfDiag(inputName = 'tmp1', outputName = 'partialPostSD', n1 = n, n2 = m, h1 = h_n, h2 = h_m)
                                    if(!predCovCurrent){
                                      remoteConstructCov(obs = FALSE, pred = TRUE, cross = FALSE, verbose = verbose)
                                      predCovCurrent <<- TRUE
                                    }
                                  } else
                                    cat("Prediction standard errors only calculated to be returned to the user, but 'ret' is FALSE.\n")

                                }                               
                                if(ret){
                                  if(se.fit){
                                    # I do the subtraction on the master because no way to extract diag of predC on slaves
                                    return(list(
                                             fit = collectVector("postMean", objPos = localProblemName, n = m, h = h_m),
                                             se.fit = sqrt(collectDiagonal("predC", objPos = localProblemName, n = m, h = h_m) -
                                             collectVector("partialPostSD", objPos = localProblemName, n = m, h = h_m))))
                                  } else{
                                    return(collectVector("postMean", objPos = localProblemName, n = m, h = h_m))
                                  }
                                } else invisible(NULL)
                              },
                              
                              
                              calcPostCov = function(returnDiag = TRUE, verbose = FALSE){
                                'Description:

Calculates the (posterior) predictive covariance (also known as the conditional variance (conditional given the observed values).


Arguments:

returnDiag: logical. Should the diagonal of the (posterior) covariance matrix (the individual conditional prediction variances for the prediction locations) be returned? 

verbose: logical. Should status information be printed?


Value:
numeric. The vector of variances, or NULL if "returnDiag" is FALSE.
'
                                
                                if(!postCovCurrent){
                                  if(!covCurrent){
                                    remoteConstructCov(obs = TRUE, pred = FALSE, cross = FALSE, verbose = verbose)
                                    covCurrent <<- TRUE
                                  }
                                  if(!cholCurrent){
                                    remoteCalcChol(matName = 'C', cholName = 'L', matPos = localProblemName, cholPos = localProblemName, n = n, h = h_n)
                                    cholCurrent <<- TRUE
                                  }
                                  if(!crossCovCurrent){
                                    remoteConstructCov(obs = FALSE, pred = FALSE, cross = TRUE, verbose = verbose)
                                    crossCovCurrent <<- TRUE
                                  }
                                  remoteForwardsolve(cholName = 'L', inputName = 'crossC', outputName = 'tmp1', cholPos = localProblemName, inputPos = localProblemName, n1 = n, n2 = m, h1 = h_n, h2 = h_m)
                                  remoteCrossProdMatSelf(inputName = 'tmp1', outputName = 'tmp2', n1 = n, n2 = m, h1 = h_n, h2 = h_m)
                                  if(!predCovCurrent){
                                    remoteConstructCov(obs = FALSE, pred = TRUE, cross = FALSE, verbose = verbose)
                                    predCovCurrent <<- TRUE
                                  }
                                  remoteCalc("predC", 'tmp2', `-`, "postC", input1Pos = localProblemName, outputPos = localProblemName)
                                  postCovCurrent <<- TRUE
                                  if(verbose)
                                    cat("Constructed posterior covariance.\n")
                                } else
                                  if(verbose)
                                    cat("Posterior covariance is current. Returning diagonal if requested.\n")
                                if(returnDiag){
                                  return(collectDiagonal("postC", objPos = localProblemName, n = m, h = h_m))
                                } else invisible(NULL)
                              },
                              
                              
                              simulateRealizations = function(r = 1, h_r = NULL, obs = FALSE, pred = FALSE, post = TRUE, verbose = FALSE){
                                'Description:
Simulates realizations in one of three ways (a) from the (prior) predictive distribution at the observation locations, (b) from the (prior) predictive distribution at the prediction locations, or (c) from the (posterior) predictive distribution at the prediction locations.


Arguments:

r: numeric. The number of realizations to simulate.

h_r: numeric. The block replication factor for the realizations; when NULL, a reasonable default value is chosen.

obs: logical. Should the realizations be simulated from the prior (unconditional on "data") predictive distribution at the observation locations?

pred: logical. Should the realizations be simulated from the prior (unconditional on "data") predictive distribution at the prediction locations?

post: logical. Should the realizations be simulated from the posterior (conditional on "data") predictive distribution at the prediction locations?

verbose: logical. Should status information be printed?


Details:
Given the set of parameter values contained in the "params" component of the krigeProblem object, simulates realizations. The standard use for this is to simulate from the posterior predictive distribution (i.e., conditional on the "data" component of the krigeProblem object at the prediction locations to characterize uncertainty in the predictions (when post = TRUE). Alternatively, one can simulate from the prior predictive distribution (i.e., not conditional on the data) at either the observation locations (when obs = TRUE) or the prediction locations (when pred = TRUE).


Value:
numeric. A matrix of realizations, with realizations as columns and prediction locations as rows.
'
                                
                                if(is.null(h_r))
                                  h_r <- calcH(r) 
                                
                                if(obs + pred + post != 1)
                                  stop("krigeProblem$simulateRealizations: Please choose either 'obs', 'pred', or 'post' as the type of realizations to simulate.")
                                if(obs){
                                  if(!meanCurrent){
                                    remoteConstructMean(obs = TRUE, pred = FALSE, verbose = verbose)
                                    meanCurrent <<- TRUE
                                  }
                                  meanVec <- collectVector('mean', n, h_n)
                                  if(!covCurrent){
                                    remoteConstructCov(obs = TRUE, pred = FALSE, verbose = verbose)
                                    covCurrent <<- TRUE
                                  }
                                  if(!cholCurrent){
                                    remoteCalcChol('C', 'L', matPos = localProblemName, cholPos = localProblemName, n = n, h = h_n)
                                    cholCurrent <<- TRUE
                                  }
                                  nVal <- n
                                  hVal <- h_n
                                  cholName <- 'L'
                                }
                                if(pred){
                                  if(!predMeanCurrent){
                                    remoteConstructMean(obs = FALSE, pred = TRUE, verbose = verbose)
                                    predMeanCurrent <<- TRUE
                                  }
                                  meanVec <- collectVector('predMean', m, h_m)
                                  if(!predCovCurrent){
                                    remoteConstructCov(obs = FALSE, pred = TRUE, verbose = verbose)
                                    predCovCurrent <<- TRUE
                                  }
                                  if(!predCholCurrent){
                                    remoteCalcChol('predC', 'predL', matPos = localProblemName, cholPos = localProblemName, n = m, h = h_m)
                                    cholCurrent <<- TRUE
                                  }
                                  nVal <- m
                                  hVal <- h_m
                                  cholName <- 'predL'
                                }
                                if(post){
                                  meanVec <- predict(ret = TRUE, verbose = verbose)
                                  if(!postCovCurrent)
                                    calcPostCov(returnDiag = FALSE, verbose = verbose)
                                  if(!postCholCurrent){
                                    remoteCalcChol('postC', 'postL', matPos = localProblemName, cholPos = localProblemName, n = m, h = h_m)
                                    postCholCurrent <<- TRUE
                                  }
                                  nVal <- m
                                  hVal <- h_m
                                  cholName <- 'postL'
                                }
                                
                                remoteConstructRnormMatrix('u', n1 = nVal, n2 = r, h1 = hVal, h2 = h_r)
                                remoteMultChol(cholName = cholName, 'u', outputName = 'result', cholPos = localProblemName, n1 = nVal, n2 = r, h1 = hVal, h2 = h_r)
                                devs <- collectRectangularMatrix('result', n1 = nVal, n2 = r, h1 = hVal, h2 = h_r)

                                return(meanVec + devs)
                              }
                              ))

#############################################################
# ReferenceClass class that holds problem on slave processes
#############################################################


distributedKrigeProblem <- setRefClass("distributedKrigeProblem",
                                  fields = list(
                                    matIndices = "ANY",
                                    vectorIndices = "ANY",
                                    predVectorIndices = "ANY",
                                    predMatIndices = "ANY",
                                    crossMatIndices = "ANY",
                                    
                                    meanFunction = "function",
                                    covFunction = "function",
                                    predMeanFunction = "function",
                                    crossCovFunction = "function",
                                    predCovFunction = "function",

                                    covFunctionCached = "ANY",
                                    crossCovFunctionCached = "ANY",
                                    predCovFunctionCached = "ANY",
                                    
                                    mean = "numeric",
                                    predMean = "numeric",
                                    postMean = "numeric",
                                    
                                    C = "numeric",
                                    crossC = "numeric",
                                    predC = "numeric",
                                    postC = "numeric",
                                    L = "numeric",
                                    predL = "numeric",
                                    postL = "numeric",

                                    data = "ANY", 
                                    params = "ANY",

                                    inputs = "list"
                                    ),
                                  
                                  methods = list(

                                    initialize = function(h_n, h_m, n, m, meanFunction = function(){}, predMeanFunction = function(){}, covFunction = function(){}, crossCovFunction = function(){}, predCovFunction = function(){}, inputs = NULL, params = NULL, ...){
                                      'Description:
Constructor that initializes the distributedKrigeProblem object on the slave processes. The object contains the metadata and distributed objects necessary for carrying out the krigeProblem methods.


Arguments:
All arguments are analogous to (and a subset of) those in the constructor for the krigeProblem object.


Value:
distributedKrigeProblem. The initialized object.
'
                                      params <<- params

                                      meanFunction <<- meanFunction
                                      predMeanFunction <<- predMeanFunction
                                      covFunction <<- covFunction
                                      crossCovFunction <<- crossCovFunction
                                      predCovFunction <<- predCovFunction
                                      
                                      covFunctionCached <<- NULL
                                      predCovFunctionCached <<- NULL
                                      crossCovFunctionCached <<- NULL
                                      
                                      inputs <<- inputs

                                      matIndices <<- localGetTriangularMatrixIndices(n, h_n)
                                      vectorIndices <<- localGetVectorIndices(n, h_n)
                                      
                                      if(m != 0){
                                        crossMatIndices <<- localGetRectangularMatrixIndices(n, m, h_n, h_m)
                                        predMatIndices <<- localGetTriangularMatrixIndices(m, h_m)
                                        predVectorIndices <<- localGetVectorIndices(m, h_m)
                                      }
                                      
                                      callSuper(...)
                                    }
                                    )
                                  )



localKrigeProblemConstructMean <- function(problemName, obs, pred){
  status <- try( {
    dProblem <- eval(as.name(problemName))
    if(pred) 
      if(.bigGP$I == .bigGP$J) 
        dProblem$predMean <- dProblem$predMeanFunction(
                               dProblem$params, dProblem$inputs, dProblem$predVectorIndices) else dProblem$predMean <- numeric(0)
    if(obs)
      if(.bigGP$I == .bigGP$J)
        dProblem$mean <- dProblem$meanFunction(dProblem$params,
                                               dProblem$inputs, dProblem$vectorIndices) else dProblem$mean <- numeric(0)
  } )
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}


localKrigeProblemConstructCov <- function(problemName, obs, pred, cross){
  status <- try( {
    dProblem <- eval(as.name(problemName))
    if(obs){
      tmp <- dProblem$covFunction(dProblem$params, dProblem$inputs, dProblem$matIndices, dProblem$covFunctionCached)
      if(is.list(tmp)){
        if(length(tmp) == 2){
          dProblem$covFunctionCached <- tmp[[2]]
          tmp <- tmp[[1]]
        } else stop("localKrigeProblemConstructCov: 'covFunction' should return a two-element list, with the covariance matrix as the first element and cached object as the second element.")
      }
      dProblem$C <- tmp
    }
    if(pred){
      tmp <- dProblem$predCovFunction(dProblem$params, dProblem$inputs, dProblem$predMatIndices, dProblem$predCovFunctionCached)
      if(is.list(tmp)){
        if(length(tmp) == 2){
          dProblem$predCovFunctionCached <- tmp[[2]]
          tmp <- tmp[[1]]
        } else stop("localKrigeProblemConstructCov: 'predCovFunction' should return a two-element list, with covariance matrix as the first element and cached object as the second element.")
      }
      dProblem$predC <- tmp
    } 
    if(cross){
      tmp <- dProblem$crossCovFunction(dProblem$params, dProblem$inputs, dProblem$crossMatIndices, dProblem$crossCovFunctionCached)
      if(is.list(tmp)){
        if(!is.null(names(tmp))){
          if(length(tmp) == 2){
            dProblem$crossCovFunctionCached <- tmp[[2]]
            tmp <- tmp[[1]]
          } else stop("localKrigeProblemConstructCov: 'crossCovFunction' should return a two-element list, with covariance matrix as the first element and cached object as the second element.")
        }
      }
      dProblem$crossC <- tmp 
    }
  } )
  if(class(status) == "try-error") invisible(status) else invisible(NULL)
}

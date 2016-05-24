####################################################
# Model based optimization based on Kriging with EGO

mbo1d <- function (model, fun, nsteps, lower, upper, parinit, isoInput, maxRuns=3, repetitions=5, tol_input=.Machine$double.eps^0.25, 
                   addInfo=TRUE, nCores=1, envir=parent.frame()) {
  
  intervalMat <- matrix(c(lower, upper), nrow=2, ncol=length(lower), byrow=TRUE)
  for(i in 1:nsteps) {
    
    # Optimize expected improvement criterion
    opt_res <- optimize1dMulti (f_input=function (x) -EImod(x=x, model=model), interval_matrix=intervalMat, maxRuns=maxRuns, 
                                repetitions=repetitions, tol_input=tol_input, x_0=parinit, addInfo=FALSE, 
                                nCores=1, envir=parent.frame(), directUse=FALSE, OptTypePar="mbo1d")
    
    # Construct new design
    new_func_val <- c(model@y, fun(opt_res$minimum))
    new_design <- rbind(model@X, opt_res$minimum)
    parinit <- new_design [which.min(new_func_val), ]
    
    # Reevaluate Kriging model
    model <- km(design=new_design, response=new_func_val, control=list(trace=FALSE), nugget.estim=TRUE, iso=isoInput)
    if(addInfo) {cat("N_step", i, "of", nsteps, "\n")}
  }
  
  # Return best solution
  Index <- which.min(new_func_val)
  result <- new_design [Index, ]
  attributes(result) <- list(funcVal=new_func_val [Index])
  return(result)
}

# Guidelines:
# Use default matern 5/2 correlation structure to avoid numerical problems during invertation of matrices
# Use number of steps and starting points proportional to dimension D, e. g. nsteps= startPoints = D * 10

# Input Arguments
# loss_func: Loss function with only a multivariate vector as input and a scalar output value!
# n_steps: Number of initial points and number of steps of EGO
# lower_bounds: Vector of lower bounds of variables
# upper_bounds: Vector of upper bounds of variables
# x_start: Starting value for optimization
# Output list with names
# par: Best parameters
# value: Function value at the best parameters

mboAll <- function (loss_func, n_steps, initDesign, lower_bounds, upper_bounds, x_start, isoInput=FALSE, addInfo=TRUE,
                    maxRuns=3, repetitions=5, tol_input=.Machine$double.eps^0.25, nCores=1, envir=parent.frame(), 
                    metaModelSelect=TRUE) {
  if(nCores<1){stop("Please specify an integer number greater or equal to one as the number of threads!")}
  # Generate design and transform 
  dimVec <- length(x_start)
  LHS_design <- maximinLHS(n=initDesign, k=dimVec)
  LHS_design <- sapply(1:dimVec, function (x) (upper_bounds[x] - lower_bounds[x]) * LHS_design [, x] + lower_bounds[x])
  if(addInfo) {cat("LHS_design", "\n")}
  
  # Evaluate loss function on design
  func_start <- vector("numeric", initDesign)
  localEnvir <- environment()
  if(nCores==1) {
    for(i in 1:initDesign) {
      func_start[i] <- loss_func(LHS_design[i, ])
      if(addInfo) {cat("Iter", i, "of", initDesign, "\n")}
    }
  }
  else {
    cl <- makeCluster(nCores)
    clusterExport(cl = cl, varlist=c("loss_func", "LHS_design", "initDesign"), envir = localEnvir)
    # Necessary to access variables in function
    tempVec <- x_start
    tempVec[seq(1, length(tempVec), 3)] <- 1
    loss_func(tempVec)[1]
    func_start <- parSapplyLB(cl=cl, X=1:initDesign, 
                              FUN=function(x) loss_func(LHS_design[x, ])[1])
    stopCluster(cl=cl)
  }
  if(addInfo) {cat("LHS_design func eval", "\n")}
  
  if(metaModelSelect) {
  # Metamodel validation and Kriging estimation
  # Trend: constant
  # covtype: "gauss", "matern5_2", "matern3_2", "exp" or "powexp"
  covtypeInput <- c("gauss", "matern5_2", "matern3_2", "exp", "powexp")
  perfMeasure <- vector("numeric", length=5)
  km_select <- vector("list", 5)
  if(nCores==1) {
    for(j in 1:5) {
      km_select[[j]] <- km(covtype=covtypeInput[j], design=LHS_design, response=func_start, control=list(trace=FALSE), nugget.estim=TRUE, iso=isoInput)
      LOOmuSigma <- leaveOneOut.km(model=km_select[[j]], type="UK", trend.reestim=TRUE)
      perfMeasure [j] <- predLogProb (predMean=LOOmuSigma$mean, predSigma=LOOmuSigma$sd^2, y=func_start, X=LHS_design)
    }
    km_base <- km_select [[which.max(perfMeasure)]]
  }
  else {
    tempFunction <- function(j) {
      km_select1 <- km(covtype=covtypeInput[j], design=LHS_design, response=func_start, control=list(trace=FALSE), nugget.estim=TRUE, iso=isoInput)
      LOOmuSigma <- leaveOneOut.km(model=km_select1, type="UK", trend.reestim=TRUE)
      perfMeasure1 <- predLogProb (predMean=LOOmuSigma$mean, predSigma=LOOmuSigma$sd^2, y=func_start, X=LHS_design)
      return(list(Perf=perfMeasure1, Model=km_select1))
    }
    cl <- makeCluster(ifelse(nCores<=5, nCores, 5))
    clusterExport(cl = cl, varlist=c("tempFunction", "func_start", "isoInput", "covtypeInput"), envir = localEnvir)
    tempResults <- parLapply(cl=cl, X=1:5, fun=tempFunction)
    stopCluster(cl=cl)
    maxInd <- which.max(sapply(1:5, function(x) tempResults[[x]]$Perf))
    km_base <- tempResults [[maxInd]]$Model
  }
  if(addInfo) {cat("Kriging model selection", "\n")}
  }
  else{
    km_base <- km(design=LHS_design, response=func_start, control=list(trace=FALSE), nugget.estim=TRUE, iso=isoInput)
  }

  # Optimize
  ego_result <- mbo1d(model=km_base, fun=loss_func, nsteps=n_steps, lower=lower_bounds, upper=upper_bounds, parinit=x_start, 
                      isoInput=isoInput, addInfo=TRUE, nCores=nCores, 
                      maxRuns=maxRuns, repetitions=repetitions, tol_input=tol_input)
  if(addInfo){cat("Mbo1D", "\n")}
  # Output
  return(list(par=c(ego_result), value=attr(ego_result, "funcVal")))
}

# Main tuning function of KDSN
tuneMboKDSN <- function (y, X, maxLevels=10, fineTuneIt=100, 
                         gammaPar=1, nStepMult=20, designMult=10, 
                         dimMax=round(sqrt(dim(X)[1])/2), addInfo=TRUE, nCores=1,
                         maxRuns=3, repetitions=5, tol_input=.Machine$double.eps^0.25) {
  
  if(nCores<1){stop("Please specify an integer number greater or equal to one as the number of threads!")}
  
#   if(nCores>1) {
#     # Initialize cluster
#     cl <- makeCluster(nCores)
#   }
  
  # Initialize parameters
  n <- dim(X)[1]
  levels <- 1
  condWhile <- TRUE
  Loss_prev <- Inf
  
  # Initialize starting vector of hyperparameters
  MaxMatEntries <- .Machine$integer.max
  dimStart <- round ((dimMax+1) / 2)
  if((dimMax*2*n) > MaxMatEntries) {
    dimMax <- floor(MaxMatEntries/n/2)
    dimStart <- round (sqrt(dimMax*2)/2)
  }
  quantEuklid <- quantile(c(dist(robustStandard(X))^2), probs=c(0, 0.5, 1))
  sigmaStart <- quantEuklid [2]
  lambdaStart <- 0
  x_start <- c(dimStart, sigmaStart, lambdaStart)
  x_new <- x_start
  x_prev <- x_start
  
  # Initialize bounds of hyperparameters
  interval_matrix_start <- matrix(c(1, dimMax, quantEuklid [1], quantEuklid [3], 0, 10), nrow=2, ncol=3)
  interval_matrix <- interval_matrix_start
  
  while(condWhile) {
    if(levels>1) {
      # Specification in loss function must match!
      # Reset starting value
      x_new <- c(x_new, x_start)
      # Compute bound matrix: rows=1 (min), 2 (max); columns=number of parameters
      interval_matrix <- cbind(interval_matrix, interval_matrix_start)
    }
    
    # Tune
    optVal <- mboAll (loss_func=function (x) lossKDSN (parOpt=x, y=y, X=X, gammaPar=gammaPar, seedW=seq(0, (levels-1), 1)), 
                      n_steps=(length(x_new)+3)*nStepMult, initDesign=(length(x_new)+3)*designMult, 
                      lower_bounds=interval_matrix [1, ], upper_bounds=interval_matrix [2, ], x_start=x_new,
                      nCores=nCores, maxRuns=maxRuns, repetitions=repetitions, tol_input=tol_input,
                      addInfo=addInfo)
    
    # Set condition
    condWhile <- optVal$value < Loss_prev
    
    # Update
    if(condWhile) {
      x_new <- optVal$par
      Loss_prev <- optVal$value
      x_prev <- x_new
      if(levels >= maxLevels) {
        break
        if(addInfo) {cat("tuneKDSN", "Level =", levels, "\n")}
      }
      levels <- levels + 1
      if(addInfo) {cat("tuneKDSN", "Level =", levels-1, "\n")}
    }
    else {
      x_new <- x_prev
      levels <- levels - 1
    }
  }
  
  # Fine tune random Fourier transformation weights
  fineTune <- vector("numeric", fineTuneIt)
  seedGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * 
                      sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), nrow=fineTuneIt, ncol=levels)
  localEnvir <- environment()
  if(nCores==1) {
    for(i in 1:fineTuneIt) {
      fineTune[i] <- lossKDSN(parOpt=x_new, y=y, X=X, gammaPar=gammaPar, seedW=seedGuess [i, ])[1]
      if(addInfo) {cat("tuneKDSN", "FineTune =", i, "\n")}
    }
  }
  else {
      cl <- makeCluster(nCores)
      clusterExport(cl = cl, varlist=c("lossKDSN", "x_new", "y", "X", "gammaPar", "seedGuess", "fineTuneIt"), 
                    envir = localEnvir)
      fineTune <- parSapply(cl=cl, X=1:fineTuneIt, 
                                FUN=function(i) lossKDSN(parOpt=x_new, y=y, X=X, gammaPar=gammaPar, seedW=seedGuess [i, ])[1])
      stopCluster(cl=cl)
      if(addInfo){cat("tuneKDSN", "FineTune done", "\n")}
  }
  minIndex <- which.min(fineTune)
  
  # Output
  # Refit best model
  lenx_new <- length(x_new)
  stopifnot((lenx_new%%3) == 0)
  levels1 <- lenx_new/3
  stopifnot(length(seedGuess[minIndex, ]) == levels1)
  Dim1 <- round(x_new[seq(1, lenx_new, 3)])
  sigma1 <- x_new[seq(2, lenx_new, 3)]
  lambda1 <- x_new[seq(3, lenx_new, 3)]
  if(fineTune[minIndex] < optVal$value){
    finalModel <- fitKDSN(y = y, X = X, levels = levels1, Dim = Dim1, 
                          sigma = sigma1, lambda = lambda1, alpha = rep(0, levels1), 
                          info = FALSE, 
                          seedW = seedGuess [minIndex, ], 
                          standX = TRUE)
  }
  else{
    finalModel <- fitKDSN(y = y, X = X, levels = levels1, Dim = Dim1, 
                          sigma = sigma1, lambda = lambda1, alpha = rep(0, levels1), 
                          info = FALSE, seedW = seq(0, (levels1-1), 1), standX = TRUE)
  }
  # Include GCV score as attribute
  attr(finalModel, which="GCV") <- fineTune[minIndex]
  return(finalModel)
}

# Main tuning function of KDSN
tuneMboCvKDSN <- function (y, X, maxLevels=10, fineTuneIt=100, 
                         nStepMult=20, designMult=10, 
                         dimMax=round(sqrt(dim(X)[1])/2), addInfo=TRUE, nCores=1,
                         maxRuns=3, repetitions=5, tol_input=.Machine$double.eps^0.25,
                         cvIndex, lossFunc=devStandard) {
  
  if(nCores<1){stop("Please specify an integer number greater or equal to one as the number of threads!")}
  
  #   if(nCores>1) {
  #     # Initialize cluster
  #     cl <- makeCluster(nCores)
  #   }
  
  # Initialize parameters
  n <- dim(X)[1]
  levels <- 1
  condWhile <- TRUE
  Loss_prev <- Inf
  
  # Initialize starting vector of hyperparameters
  MaxMatEntries <- .Machine$integer.max
  dimStart <- round ((dimMax+1) / 2)
  if((dimMax*2*n) > MaxMatEntries) {
    dimMax <- floor(MaxMatEntries/n/2)
    dimStart <- round (sqrt(dimMax*2)/2)
  }
  quantEuklid <- quantile(c(dist(robustStandard(X))^2), probs=c(0, 0.5, 1))
  sigmaStart <- quantEuklid [2]
  lambdaStart <- 0
  x_start <- c(dimStart, sigmaStart, lambdaStart)
  x_new <- x_start
  x_prev <- x_start
  
  # Initialize bounds of hyperparameters
  interval_matrix_start <- matrix(c(1, dimMax, quantEuklid [1], quantEuklid [3], 0, 10), nrow=2, ncol=3)
  interval_matrix <- interval_matrix_start
  
  while(condWhile) {
    if(levels>1) {
      # Specification in loss function must match!
      # Reset starting value
      x_new <- c(x_new, x_start)
      # Compute bound matrix: rows=1 (min), 2 (max); columns=number of parameters
      interval_matrix <- cbind(interval_matrix, interval_matrix_start)
    }
    
    # Tune
    optVal <- mboAll (loss_func=function (x) lossCvKDSN (parOpt=x, y=y, X=X, cvIndex=cvIndex,  
                                                       seedW=seq(0, (levels-1), 1), lossFunc=lossFunc), 
                      n_steps=(length(x_new)+3)*nStepMult, initDesign=(length(x_new)+3)*designMult, 
                      lower_bounds=interval_matrix [1, ], upper_bounds=interval_matrix [2, ], x_start=x_new,
                      nCores=nCores, maxRuns=maxRuns, repetitions=repetitions, tol_input=tol_input,
                      addInfo=addInfo)
    
    # Set condition
    condWhile <- optVal$value < Loss_prev
    
    # Update
    if(condWhile) {
      x_new <- optVal$par
      Loss_prev <- optVal$value
      x_prev <- x_new
      if(levels >= maxLevels) {
        break
        if(addInfo) {cat("tuneKDSN", "Level =", levels, "\n")}
      }
      levels <- levels + 1
      if(addInfo) {cat("tuneKDSN", "Level =", levels-1, "\n")}
    }
    else {
      x_new <- x_prev
      levels <- levels - 1
    }
  }
  
  # Fine tune random Fourier transformation weights
  fineTune <- vector("numeric", fineTuneIt)
  seedGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * 
                        sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), nrow=fineTuneIt, ncol=levels)
  localEnvir <- environment()
  if(nCores==1) {
    for(i in 1:fineTuneIt) {
      fineTune[i] <- lossCvKDSN(parOpt=x_new, y=y, X=X, seedW=seedGuess [i, ], 
                                cvIndex=cvIndex, lossFunc=lossFunc)[1]
      if(addInfo) {cat("tuneKDSN", "FineTune =", i, "\n")}
    }
  }
  else {
    cl <- makeCluster(nCores)
    clusterExport(cl = cl, varlist=c("lossCvKDSN", "x_new", "y", "X", "seedGuess", "fineTuneIt"), 
                  envir = localEnvir)
    fineTune <- parSapply(cl=cl, X=1:fineTuneIt, 
                          FUN=function(i) lossCvKDSN(parOpt=x_new, y=y, X=X, seedW=seedGuess [i, ],
                                                     cvIndex=cvIndex, lossFunc=lossFunc)[1])
    stopCluster(cl=cl)
    if(addInfo){cat("tuneKDSN", "FineTune done", "\n")}
  }
  minIndex <- which.min(fineTune)
  
  # Output
  # Refit best model
  lenx_new <- length(x_new)
  stopifnot((lenx_new%%3) == 0)
  levels1 <- lenx_new/3
  stopifnot(length(seedGuess[minIndex, ]) == levels1)
  Dim1 <- round(x_new[seq(1, lenx_new, 3)])
  sigma1 <- x_new[seq(2, lenx_new, 3)]
  lambda1 <- x_new[seq(3, lenx_new, 3)]
  if(fineTune[minIndex] < optVal$value){
    finalModel <- fitKDSN(y = y, X = X, levels = levels1, Dim = Dim1, 
                          sigma = sigma1, lambda = lambda1, alpha = rep(0, levels1), 
                          info = FALSE, 
                          seedW = seedGuess [minIndex, ], 
                          standX = TRUE)
  }
  else{
    finalModel <- fitKDSN(y = y, X = X, levels = levels1, Dim = Dim1, 
                          sigma = sigma1, lambda = lambda1, alpha = rep(0, levels1), 
                          info = FALSE, seedW = seq(0, (levels1-1), 1), standX = TRUE)
  }
  # Include GCV score as attribute
  attr(finalModel, which="GCV") <- fineTune[minIndex]
  return(finalModel)
}

# Main tuning function of KDSN
tuneMboLevelKDSN <- function (y, X, levels=10, fineTuneIt=100, gammaPar=1, 
                              nStepMult=20, designMult=10, dimMax=round(sqrt(dim(X)[1])/2), addInfo=TRUE, nCores=1,
                              maxRuns=3, repetitions=5, tol_input=.Machine$double.eps^0.25) {
  
  if(nCores<1){stop("Please specify an integer number greater or equal to one as the number of threads!")}
#   if(nCores>1) {
#     # Initialize cluster
#     cl <- makeCluster(nCores)
#   }
   
  # Initialize parameters
  n <- dim(X)[1]

  # Initialize starting vector of hyperparameters
  MaxMatEntries <- .Machine$integer.max
  dimStart <- round ((dimMax+1) / 2)
  if((dimMax*2*n) > MaxMatEntries) {
    dimMax <- floor(MaxMatEntries/n/2)
    dimStart <- round (sqrt(dimMax*2)/2)
  }
  quantEuklid <- quantile(c(dist(robustStandard(X))^2), probs=c(0, 0.5, 1))
  sigmaStart <- quantEuklid [2]
  lambdaStart <- 0
  x_start <- c(dimStart, sigmaStart, lambdaStart)
  x_new <- rep(x_start, levels)

  # Initialize bounds of hyperparameters
  interval_matrix_start <- matrix(c(1, dimMax, quantEuklid [1], quantEuklid [3], 0, 10), nrow=2, ncol=3)
  interval_matrix <- interval_matrix_start [, rep(1:dim(interval_matrix_start)[2], levels)]
  
  # Tune
  optVal <- mboAll (loss_func=function (x) lossKDSN (parOpt=x, y=y, X=X, gammaPar=gammaPar, seedW=seq(0, (levels-1), 1)), 
                           n_steps=(length(x_new)+3)*nStepMult, initDesign=(length(x_new)+3)*designMult, 
                    lower_bounds=interval_matrix [1, ], upper_bounds=interval_matrix [2, ], x_start=x_new,
                    nCores=nCores, addInfo=addInfo, 
                    maxRuns=maxRuns, repetitions=repetitions, tol_input=tol_input)
  x_new <- optVal$par

  # Fine tune random fourier transformation weights
  # Reproduceability is ensured with seed generation
  fineTune <- vector("numeric", fineTuneIt)
  seedGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), 
                      nrow=fineTuneIt, ncol=levels)
  localEnvir <- environment()
  if(nCores==1) {
    for(i in 1:fineTuneIt) {
      fineTune[i] <- lossKDSN(parOpt=x_new, y=y, X=X, gammaPar=gammaPar, seedW=seedGuess [i, ])[1]
      if(addInfo) {cat("tuneKDSN", "FineTune =", i, "\n")}
    }
  }
  else {
    cl <- makeCluster(nCores)
    clusterExport(cl = cl, varlist=c("lossKDSN", "x_new", "y", "X", "gammaPar", "seedGuess", "fineTuneIt"), 
                  envir = localEnvir)
    fineTune <- parSapply(cl=cl, X=1:fineTuneIt, 
                            FUN=function(i) lossKDSN(parOpt=x_new, y=y, X=X, gammaPar=gammaPar, seedW=seedGuess [i, ])[1])
    stopCluster(cl=cl)
    if(addInfo) {cat("FineTuning", "\n")}
  }
  minIndex <- which.min(fineTune)
  
  # Output
  # Refit best model
  lenx_new <- length(x_new)
  stopifnot((lenx_new%%3) == 0)
  levels1 <- lenx_new/3
  stopifnot(length(seedGuess[minIndex, ]) == levels1)
  Dim1 <- round(x_new[seq(1, lenx_new, 3)])
  sigma1 <- x_new[seq(2, lenx_new, 3)]
  lambda1 <- x_new[seq(3, lenx_new, 3)]
  if(fineTune[minIndex] < optVal$value){
    finalModel <- fitKDSN(y = y, X = X, levels = levels1, Dim = Dim1, 
                          sigma = sigma1, lambda = lambda1, alpha = rep(0, levels1), 
                          info = FALSE, 
                          seedW = seedGuess [minIndex, ], 
                          standX = TRUE)
    # Include GCV score as attribute
    attr(finalModel, which="GCV") <- fineTune[minIndex]
  }
  else{
    finalModel <- fitKDSN(y = y, X = X, levels = levels1, Dim = Dim1, 
                          sigma = sigma1, lambda = lambda1, alpha = rep(0, levels1), 
                          info = FALSE, 
                          seedW = seq(0, (levels1-1), 1), standX = TRUE)
    # Include GCV score as attribute
    attr(finalModel, which="GCV") <- optVal$value
  }
  return(finalModel)
}

# Main tuning function of KDSN
tuneMboLevelCvKDSN <- function (y, X, levels=10, fineTuneIt=100, nStepMult=20, designMult=10, 
                                dimMax=round(sqrt(dim(X)[1])/2), addInfo=TRUE, nCores=1,
                              maxRuns=3, repetitions=5, tol_input=.Machine$double.eps^0.25, 
                              cvIndex, lossFunc=devStandard) {
  
  if(nCores<1){stop("Please specify an integer number greater or equal to one as the number of threads!")}
  #   if(nCores>1) {
  #     # Initialize cluster
  #     cl <- makeCluster(nCores)
  #   }
  
  # Initialize parameters
  n <- dim(X)[1]
  
  # Initialize starting vector of hyperparameters
  MaxMatEntries <- .Machine$integer.max
  dimStart <- round ((dimMax+1) / 2)
  if((dimMax*2*n) > MaxMatEntries) {
    dimMax <- floor(MaxMatEntries/n/2)
    dimStart <- round (sqrt(dimMax*2)/2)
  }
  quantEuklid <- quantile(c(dist(robustStandard(X))^2), probs=c(0, 0.5, 1))
  sigmaStart <- quantEuklid [2]
  lambdaStart <- 0
  x_start <- c(dimStart, sigmaStart, lambdaStart)
  x_new <- rep(x_start, levels)
  
  # Initialize bounds of hyperparameters
  interval_matrix_start <- matrix(c(1, dimMax, quantEuklid [1], quantEuklid [3], 0, 10), nrow=2, ncol=3)
  interval_matrix <- interval_matrix_start [, rep(1:dim(interval_matrix_start)[2], levels)]
  
  # Tune
  optVal <- mboAll (loss_func=function (x) lossCvKDSN (parOpt=x, y=y, X=X, cvIndex=cvIndex,  
                                                     seedW=seq(0, (levels-1), 1), lossFunc=lossFunc), 
                    n_steps=(length(x_new)+3)*nStepMult, initDesign=(length(x_new)+3)*designMult, 
                    lower_bounds=interval_matrix [1, ], upper_bounds=interval_matrix [2, ], x_start=x_new,
                    nCores=nCores, addInfo=addInfo, 
                    maxRuns=maxRuns, repetitions=repetitions, tol_input=tol_input)
  x_new <- optVal$par
  
  # Fine tune random fourier transformation weights
  # Reproduceability is ensured with seed generation
  fineTune <- vector("numeric", fineTuneIt)
  seedGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), 
                      nrow=fineTuneIt, ncol=levels)
  localEnvir <- environment()
  if(nCores==1) {
    for(i in 1:fineTuneIt) {
      fineTune[i] <- lossCvKDSN(parOpt=x_new, y=y, X=X, cvIndex=cvIndex, 
                                seedW=seedGuess [i, ], lossFunc=lossFunc)[1]
      if(addInfo) {cat("tuneKDSN", "FineTune =", i, "\n")}
    }
  }
  else {
    cl <- makeCluster(nCores)
    clusterExport(cl = cl, varlist=c("lossCvKDSN", "x_new", "y", "X", "seedGuess", "fineTuneIt"), 
                  envir = localEnvir)
    fineTune <- parSapply(cl=cl, X=1:fineTuneIt, 
                          FUN=function(i) lossCvKDSN(parOpt=x_new, y=y, X=X, cvIndex=cvIndex,  
                                                     seedW=seedGuess [i, ], lossFunc=lossFunc)[1])
    stopCluster(cl=cl)
    if(addInfo) {cat("FineTuning", "\n")}
  }
  minIndex <- which.min(fineTune)
  
  # Output
  # Refit best model
  lenx_new <- length(x_new)
  stopifnot((lenx_new%%3) == 0)
  levels1 <- lenx_new/3
  stopifnot(length(seedGuess[minIndex, ]) == levels1)
  Dim1 <- round(x_new[seq(1, lenx_new, 3)])
  sigma1 <- x_new[seq(2, lenx_new, 3)]
  lambda1 <- x_new[seq(3, lenx_new, 3)]
  if(fineTune[minIndex] < optVal$value){
    finalModel <- fitKDSN(y = y, X = X, levels = levels1, Dim = Dim1, 
                          sigma = sigma1, lambda = lambda1, alpha = rep(0, levels1), 
                          info = FALSE, 
                          seedW = seedGuess [minIndex, ], 
                          standX = TRUE)
    attr(finalModel, which="Loss") <- fineTune[minIndex]
  }
  else{
    finalModel <- fitKDSN(y = y, X = X, levels = levels1, Dim = Dim1, 
                          sigma = sigma1, lambda = lambda1, alpha = rep(0, levels1), 
                          info = FALSE, 
                          seedW = seq(0, (levels1-1), 1), standX = TRUE)
    attr(finalModel, which="Loss") <- optVal$value
  }
  # Include Loss score, loss function and cvIndex as attributes
  attr(finalModel, which="LossFunc") <- lossFunc
  attr(finalModel, which="cvIndex") <- cvIndex
  return(finalModel)
}

########################################
# Reformulate EI without numerical check
# Check is flawed and gives error in some cases

EImod <- function (x, model, plugin = NULL, type = "UK", minimization = TRUE, 
                   envir = NULL) 
{
  if (is.null(plugin)) {
    if (minimization) {
      plugin <- min(model@y)
    }
    else {
      plugin <- -max(model@y)
    }
  }
  m <- plugin
  d <- length(x)
  if (d != model@d) {
    stop("x does not have the right size")
  }
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  predx <- predict(object = model, newdata = newdata, type = type, 
                   checkNames = FALSE)
  kriging.mean <- predx$mean
  if (!minimization) {
    kriging.mean <- -kriging.mean
  }
  kriging.sd <- predx$sd
  xcr <- (m - kriging.mean)/kriging.sd
  #######################
  # Code in function EI()
  #   if (kriging.sd/sqrt(model@covariance@sd2) < 1e-06) {
  #     res <- 0
  #     xcr <- xcr.prob <- xcr.dens <- NULL
  #   }
  # else{
  #  xcr.prob <- pnorm(xcr)
  #  xcr.dens <- dnorm(xcr)
  #  res <- (m - kriging.mean) * xcr.prob + kriging.sd * xcr.dens
  # }
  
  xcr.prob <- pnorm(xcr)
  xcr.dens <- dnorm(xcr)
  res <- (m - kriging.mean) * xcr.prob + kriging.sd * xcr.dens
  
  # Numerical check of kriging.mean, kriging.sd inputs in EI formula
  CheckCondition <- is.infinite(res) | 
                    is.null(res) | 
                    is.nan(res) |
                    is.na(res)
  if(any(CheckCondition)) {
    res[CheckCondition] <- 0
  }
  
  if (!is.null(envir)) {
    assign("xcr", xcr, envir = envir)
    assign("xcr.prob", xcr.prob, envir = envir)
    assign("xcr.dens", xcr.dens, envir = envir)
    assign("kriging.sd", kriging.sd, envir = envir)
    assign("c", predx$c, envir = envir)
    assign("Tinv.c", predx$Tinv.c, envir = envir)
  }
  return(res)
}

###############################################
# Grid search over subset of levels
# In each level MBO algorithm will be performed

tuneMboLevelGridKDSN <- function(y, X, levelSet, fineTuneIt=100, gammaPar=1, 
                                 nStepMult=20, designMult=10, dimMax=round(sqrt(dim(X)[1])/2), addInfo=TRUE, 
                                 nCoresInner=1, nCoresOuter=1, maxRuns=3, repetitions=5, 
                                 tol_input=.Machine$double.eps^0.25) {
  # Check parallel arguments
  if(round(nCoresInner)!=nCoresInner & round(nCoresOuter)!=nCoresOuter) {
    stop("Please specify integer numbers in nCoresInner, nCoresOuter!")}
  if(nCoresInner>1 & nCoresOuter>1) {stop("Only one parallelisation technique is allowed. 
                                          Please specify either nCoresInner or nCoresOuter 
                                          with integer numbers greater 1!")}
  
  # Compute KDSN with MBO tuning
  localEnvir <- environment()
  if(nCoresOuter==1) {
    resGridLevel <- vector("list", length(levelSet))
    for(l in levelSet) {
      resGridLevel[[l]] <- tuneMboLevelKDSN(y=y, X=X, levels=levelSet[l], fineTuneIt=fineTuneIt, gammaPar=gammaPar, 
                        nStepMult=nStepMult, designMult=designMult, dimMax=dimMax, addInfo=addInfo, 
                        nCores=nCoresInner, maxRuns=maxRuns, repetitions=repetitions, tol_input=tol_input)
      cat("Level", levelSet[l], "\n")
    }
  }
  else{
    cl <- makeCluster(nCoresOuter)
    clusterExport(cl = cl, varlist=c(ls(), "tuneMboLevelKDSN"), envir = localEnvir)
    resGridLevel <- parLapplyLB(cl=cl, X=levelSet, 
                                fun=function(x) tuneMboLevelKDSN(y=y, X=X, levels=levelSet[x], fineTuneIt=fineTuneIt, gammaPar=gammaPar, 
                                                nStepMult=nStepMult, designMult=designMult, dimMax=dimMax, addInfo=addInfo, 
                                                nCores=nCoresInner, maxRuns=maxRuns, repetitions=repetitions, tol_input=tol_input))
    stopCluster(cl=cl)
  }
  
  # Output tuned KDSN with smallest GCV
  gridLevelScores <- sapply(1:length(resGridLevel), function(x) attr(resGridLevel[[x]], "GCV"))
  return(resGridLevel[[which.min(gridLevelScores)]])
}

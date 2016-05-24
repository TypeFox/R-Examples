# General sequential one dimensional optimization of f(x), x \in R^d, f(x) \in R
optimize1dMulti <- function (f_input, interval_matrix, maxRuns=3, repetitions=5, 
                             tol_input=.Machine$double.eps^0.25, x_0=NULL, addInfo=TRUE,
                             nCores=1, envir=parent.frame(), directUse=TRUE, OptTypePar="") {
  
  # Check if first argument of function is x
  stopifnot(formalArgs (f_input) [1]=="x")

  # Rerun optimization with different starting values
  Results <- vector("list", repetitions)
  dimension <- dim(interval_matrix) [2]
  if(nCores==1) {
    for(j in 1:repetitions) {
      
      if(j > 1 | is.null(x_0)) {
        # Set initial starting value: Random vector x nid ~ U (a_i, b_i)
        x_0 <- sapply(1:dimension, function (x) runif(n=1, min=interval_matrix [1, x], max=interval_matrix [2, x]))
      }
      x_0_alt <- x_0
      liste_vektoren <- vector("list", dimension)
      
      if(dimension==1) {
        liste_vektoren [[1]] <- expression(x)
      }
      
      if(dimension==2) {
        liste_vektoren [[1]] <- expression(c(x, x_0 [2]))
        liste_vektoren [[2]] <- expression(c(x_0 [1], x))
      }
      
      if(dimension>=3) {
        liste_vektoren [[1]] <- expression(c(x, x_0 [2:dimension]))
        liste_vektoren [[dimension]] <- expression(c(x_0 [1:(dimension-1)], x))
        for(i in 1:(dimension-2)) {
          liste_vektoren [[i+1]] <- substitute(c(x_0 [1:i], x, x_0 [(2+i):dimension]),list(i=i))
        }
      }
      
      # Univariate optimization over one variable given all other variables
      i <- 1
      whileCondition <- TRUE
      stepRuns <- 0
      f_input_x_0 <- f_input(x_0)
      f_input_x_0_alt <- f_input_x_0
      x_0_alt <- x_0
      
      while (whileCondition) {
        # Conditional optimization
        for(i in 1:dimension) {
          erg <- optimize(f=function (x) f_input( eval( liste_vektoren [[i]] ) ), interval=interval_matrix [, i], tol=tol_input)
          x_0 [i] <- erg$minimum
          if(addInfo) {
            cat("optimize1dMulti", "parameter", i, "\n")
          }
        }
        
        # Condition for stopping
        stepRuns <- stepRuns + 1
        f_input_x_0 <- erg$objective
        if(stepRuns == maxRuns) {
          whileCondition <- FALSE
        }
        else {
          if(abs(f_input_x_0_alt) != 0 & sum(abs(x_0_alt)) != 0) {
            whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) / abs(f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) / sum(abs(x_0_alt)) >= tol_input)
          }
          else {
            if(abs(f_input_x_0_alt) == 0 & sum(abs(x_0_alt)) != 0) {
              whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) / sum(abs(x_0_alt)) >= tol_input)
            }
            if(abs(f_input_x_0_alt) != 0 & sum(abs(x_0_alt)) == 0) {
              whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) / abs(f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) >= tol_input)
            }
            if(abs(f_input_x_0_alt) == 0 & sum(abs(x_0_alt)) == 0) {
              whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) >= tol_input)
            }
          }
        }
        f_input_x_0_alt <- f_input_x_0
        x_0_alt <- x_0
        if(addInfo) {
          cat("optimize1dMulti", "run", stepRuns, "\n")
        }
      }
      Results [[j]] <- list(minimum=x_0, objective=f_input_x_0)
      if(addInfo) {
        cat("optimize1dMulti", "repetition", j, "\n")
      }
    }
  }
  else{

    # Help function for parallelisation
    tempFunc <- function(j) {
      
      if(j > 1 | is.null(x_0)) {
        # Set initial starting value: Random vector x nid ~ U (a_i, b_i)
        x_0 <- sapply(1:dimension, function (x) runif(n=1, min=interval_matrix [1, x], max=interval_matrix [2, x]))
      }
      x_0_alt <- x_0
      liste_vektoren <- vector("list", dimension)
      
      if(dimension==1) {
        liste_vektoren [[1]] <- expression(x)
      }
      
      if(dimension==2) {
        liste_vektoren [[1]] <- expression(c(x, x_0 [2]))
        liste_vektoren [[2]] <- expression(c(x_0 [1], x))
      }
      
      if(dimension>=3) {
        liste_vektoren [[1]] <- expression(c(x, x_0 [2:dimension]))
        liste_vektoren [[dimension]] <- expression(c(x_0 [1:(dimension-1)], x))
        for(i in 1:(dimension-2)) {
          liste_vektoren [[i+1]] <- substitute(c(x_0 [1:i], x, x_0 [(2+i):dimension]),list(i=i))
        }
      }

      # Univariate optimization of one variable given all other variables
      i <- 1
      whileCondition <- TRUE
      stepRuns <- 0
      f_input_x_0 <- f_input(x_0)
      f_input_x_0_alt <- f_input_x_0
      x_0_alt <- x_0
      
      while (whileCondition) {
        # Conditional optimization
        for(i in 1:dimension) {
          erg <- optimize(f=function (x) f_input( eval( liste_vektoren [[i]] ) ), interval=interval_matrix [, i], tol=tol_input)
          x_0 [i] <- erg$minimum
          if(addInfo) {
            cat("optimize1dMulti", "parameter", i, "\n")
          }
        }
        
        # Condition for stopping
        stepRuns <- stepRuns + 1
        f_input_x_0 <- erg$objective
        if(stepRuns == maxRuns) {
          whileCondition <- FALSE
        }
        else {
          if(abs(f_input_x_0_alt) != 0 & sum(abs(x_0_alt)) != 0) {
            whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) / abs(f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) / sum(abs(x_0_alt)) >= tol_input)
          }
          else {
            if(abs(f_input_x_0_alt) == 0 & sum(abs(x_0_alt)) != 0) {
              whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) / sum(abs(x_0_alt)) >= tol_input)
            }
            if(abs(f_input_x_0_alt) != 0 & sum(abs(x_0_alt)) == 0) {
              whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) / abs(f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) >= tol_input)
            }
            if(abs(f_input_x_0_alt) == 0 & sum(abs(x_0_alt)) == 0) {
              whileCondition <- (abs(f_input_x_0 - f_input_x_0_alt) >= tol_input) & (sum(abs(x_0 - x_0_alt)) >= tol_input)
            }
          }
        }
        f_input_x_0_alt <- f_input_x_0
        x_0_alt <- x_0
        if(addInfo) {
          cat("optimize1dMulti", "run", stepRuns, "\n")
        }
      }
      Results1 <- list(minimum=x_0, objective=f_input_x_0)
      if(addInfo) {
        cat("optimize1dMulti", "repetition", j, "\n")
      }
      return(Results1)
    }
    
    if(!directUse){
      
      # Subroutine to direct minimization
      if(OptTypePar=="tuneKDSN" ||  OptTypePar=="tuneLevelKDSN") {
        localEnvir <- environment()
        clusterExport(cl=envir$cl, varlist=ls(), envir=localEnvir)
        clusterExport(cl = envir$cl, varlist = objects(envir), envir=envir)
        Results <- parLapply(cl = envir$cl, X=1:repetitions, fun=tempFunc)
      }
      
      # Subroutine to Mbo minimization
      if(OptTypePar=="mbo1d") {
        localEnvir <- environment()
        clusterExport(cl=envir$envir$cl, varlist=ls(), envir=localEnvir)
        Results <- parLapply(cl = envir$envir$cl, X=1:repetitions, fun=tempFunc)
      }
    }
    else{
      cl <- makeCluster(nCores)
      localEnvir <- environment()
      clusterExport(cl=cl, varlist=ls(all.names=TRUE), envir=localEnvir)
      Results <- parLapply(cl = cl, X=1:repetitions, fun=tempFunc)
      stopCluster(cl=cl)
    }
  }
    
  # Choose best iteration among the random starting values
  Index <- which.min(sapply(1:repetitions, function (x) Results [[x]]$objective))
  Output <- list(minimum=Results [[Index]]$minimum, objective=Results [[Index]]$objective)
  return(Output)
}

# Main tuning function of KDSN
tuneKDSN <- function (y, X, maxRuns=3, repetitions=5, maxLevels=10, gammaPar=1, 
                      fineTuneIt=100, tol_input=.Machine$double.eps^0.25, addInfo=TRUE, dimMax=round(sqrt(dim(X)[1])/2),
                      nCores=1) {
  
  if(nCores>1) {
    # Initialize cluster
    cl <- makeCluster(nCores)
  }
  
  # Initialize parameters
  n <- dim(X) [1]
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
    optVal <- optimize1dMulti (f_input=function (x) lossKDSN (parOpt=x, 
                               y=y, X=X, gammaPar=gammaPar, seedW=seq(0, (levels-1), 1)), 
                               interval_matrix=interval_matrix, maxRuns=maxRuns, repetitions=repetitions, x_0=x_new, 
                               tol_input=tol_input, addInfo=addInfo, directUse=FALSE, nCores=nCores, OptTypePar="tuneKDSN")
    
    # Set condition
    condWhile <- optVal$objective < Loss_prev
    
    # Update
    if(condWhile) {
      x_new <- optVal$minimum
      Loss_prev <- optVal$objective
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
  
  # Fine tune random fourier transformation weights
  # Reproduceability is ensured with seed generation
  fineTune <- vector("numeric", fineTuneIt)
  seedGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * 
                        sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), nrow=fineTuneIt, ncol=levels)
  if(nCores==1) {
    for(i in 1:fineTuneIt) {
      fineTune[i] <- lossKDSN(parOpt=x_new, y=y, X=X, gammaPar=gammaPar, seedW=seedGuess [i, ])[1]
      if(addInfo) {cat("tuneKDSN", "FineTune =", i, "\n")}
    }
  }
  else {
    localEnvir <- environment()
    clusterExport(cl = cl, varlist=c("lossKDSN", "x_new", "y", "X", "gammaPar", "seedGuess", "fineTuneIt"), 
                  envir = localEnvir)
    fineTune <- parSapply(cl=cl, X=1:fineTuneIt, 
                          FUN=function(i) lossKDSN(parOpt=x_new, y=y, X=X, gammaPar=gammaPar, seedW=seedGuess [i, ])[1])
    stopCluster(cl=cl)
    if(addInfo) {cat("tuneKDSN", "FineTune done", "\n")}
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
  if(fineTune[minIndex] < optVal$objective){
    finalModel <- fitKDSN(y = y, X = X, levels = levels1, Dim = Dim1, 
                          sigma = sigma1, lambda = lambda1, alpha = rep(0, levels1), 
                          info = FALSE, 
                          seedW = seedGuess [minIndex, ], standX = TRUE)
    # Include GCV score as attribute
    attr(finalModel, which="GCV") <- fineTune[minIndex]
  }
  else{
    finalModel <- fitKDSN(y = y, X = X, levels = levels1, Dim = Dim1, 
                          sigma = sigma1, lambda = lambda1, alpha = rep(0, levels1), 
                          info = FALSE, 
                          seedW = seq(0, (levels1-1), 1), standX = TRUE)
    # Include GCV score as attribute
    attr(finalModel, which="GCV") <- optVal$objective
  }
  return(finalModel)
}

# Main tuning function of KDSN
tuneLevelKDSN <- function (y, X, maxRuns=3, repetitions=5, levels=10, gammaPar=1, 
                           fineTuneIt=100, tol_input=.Machine$double.eps^0.25, addInfo=TRUE, dimMax=round(sqrt(dim(X)[1])/2),
                           nCores=1) {
  
  if(nCores>1) {
    # Initialize cluster
    cl <- makeCluster(nCores)
  }
  
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
    optVal <- optimize1dMulti (f_input=function (x) lossKDSN (parOpt=x, 
                               y=y, X=X, gammaPar=gammaPar, seedW=seq(0, (levels-1), 1)), 
                               interval_matrix=interval_matrix, maxRuns=maxRuns, repetitions=repetitions, x_0=x_new, 
                               tol_input=tol_input, addInfo=addInfo, directUse=FALSE, 
                               OptTypePar="tuneLevelKDSN", nCores=nCores)
    x_new <- optVal$minimum
    if(addInfo){cat("tuneKDSN", "Optimize 1D done", "\n")}
    
    # Fine tune random fourier transformation weights
    # Reproduceability is ensured with seed generation
    fineTune <- vector("numeric", fineTuneIt)
    seedGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), nrow=fineTuneIt, ncol=levels)
    if(nCores==1) {
      for(i in 1:fineTuneIt) {
        fineTune[i] <- lossKDSN(parOpt=x_new, y=y, X=X, gammaPar=gammaPar, seedW=seedGuess [i, ])[1]
        if(addInfo) {cat("tuneLevelKDSN", "FineTune =", i, "\n")}
      }
    }
    else {
      localEnvir <- environment()
      clusterExport(cl = cl, varlist=c("lossKDSN", "x_new", "y", "X", "gammaPar", "seedGuess", "fineTuneIt"), 
                    envir = localEnvir)
      fineTune <- parSapply(cl=cl, X=1:fineTuneIt, 
                            FUN=function(i) lossKDSN(parOpt=x_new, y=y, X=X, gammaPar=gammaPar, seedW=seedGuess [i, ])[1])
      stopCluster(cl=cl)
      if(addInfo){cat("tuneLevelKDSN", "FineTune done", "\n")}
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
    if(fineTune[minIndex] < optVal$objective){
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
      attr(finalModel, which="GCV") <- optVal$objective
    }
  return(finalModel)
}

###############################################
# Grid search over subset of levels
# In each level MBO algorithm will be performed

tuneLevelGridKDSN <- function(y, X, maxRuns=3, repetitions=5, levelSet, gammaPar=1, 
                              fineTuneIt=100, tol_input=.Machine$double.eps^0.25, addInfo=TRUE, 
                              dimMax=round(sqrt(dim(X)[1])/2), nCoresInner=1, nCoresOuter=1) {
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
      resGridLevel[[l]] <- tuneLevelKDSN(y=y, X=X, levels=levelSet[l], fineTuneIt=fineTuneIt, gammaPar=gammaPar, 
                                        dimMax=dimMax, addInfo=addInfo, nCores=nCoresInner, maxRuns=maxRuns, 
                                        repetitions=repetitions, tol_input=tol_input)
      cat("Level", levelSet[l], "\n")
    }
  }
  else{
    cl <- makeCluster(nCoresOuter)
    clusterExport(cl = cl, varlist=c(ls(), "tuneLevelKDSN"), envir = localEnvir)
    resGridLevel <- parLapplyLB(cl=cl, X=levelSet, 
                                fun=function(x) tuneLevelKDSN(y=y, X=X, levels=levelSet[x], fineTuneIt=fineTuneIt, gammaPar=gammaPar, 
                                                              dimMax=dimMax, addInfo=addInfo, nCores=nCoresInner, maxRuns=maxRuns, 
                                                              repetitions=repetitions, tol_input=tol_input))
    stopCluster(cl=cl)
  }
  
  # Output tuned KDSN with smallest GCV
  gridLevelScores <- sapply(1:length(resGridLevel), function(x) attr(resGridLevel[[x]], "GCV"))
  return(resGridLevel[[which.min(gridLevelScores)]])
}

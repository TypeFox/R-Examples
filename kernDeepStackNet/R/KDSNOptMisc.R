################################
# Finetuning of given KDSN model

fineTuneKDSN <- function(estKDSN, y, X, gammaPar=1, fineTuneIt=100, info=TRUE, 
                         seedInit=NULL, ncpus=1) {
  # Input adjustments
  if(ncpus<1 | floor(ncpus)!=ncpus) {stop("Please supply a positive number of available cpus greater 
                    or equal to 1 (integer scalar)")}
  fineTune <- vector("numeric", fineTuneIt)
  levels <- estKDSN$Input$levels
  x_new <- c(matrix(c(estKDSN$Input$Dim, estKDSN$Input$sigma, estKDSN$Input$lambda), byrow=TRUE, nrow=3))
  
  # Reproduceability is ensured with seed generation
  set.seed(seedInit)
  seedGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), nrow=fineTuneIt, ncol=levels)
  envir1 <- environment()
  
  if(ncpus==1) {
    # Fine tune random Fourier transformation weights
    for(i in 1:fineTuneIt) {
      fineTune [i] <- lossKDSN (parOpt=x_new, y=y, X=X, gammaPar=gammaPar, seedW=seedGuess [i, ])
      if(info) {cat("tuneKDSN", "FineTune =", i, "\n")}
    }
  }
  else {
    tempFunc <- function(i) {
      fineTuneResult <- lossKDSN (parOpt=x_new, y=y, X=X, gammaPar=gammaPar, seedW=seedGuess [i, ])
      return(fineTuneResult)
    }
    
    # Parallel execution with parallel package
    cl <- makeCluster(ncpus)
    clusterExport(cl=cl, varlist=c("seedGuess", "x_new", "y", "X", "gammaPar"), 
                  envir=envir1)
    clusterEvalQ(cl=cl, library(kernDeepStackNet))
    fineTune <- parSapply(cl=cl, X=1:fineTuneIt, FUN=tempFunc)
    stopCluster(cl)
  }
  minIndex <- which.min(fineTune)

  # Evaluate loss of given model
  givenLoss <- c(lossKDSNgivenModel (KDSNfit=estKDSN, y=y, X=X, gammaPar=gammaPar))
  
  # Refit best model and output if it is better than the given model
  lenx_new <- length(x_new)
  stopifnot((lenx_new%%3) == 0)
  levels1 <- lenx_new/3
  stopifnot(length(seedGuess[minIndex, ]) == levels1)
  Dim1 <- round(x_new[seq(1, lenx_new, 3)])
  sigma1 <- x_new[seq(2, lenx_new, 3)]
  lambda1 <- x_new[seq(3, lenx_new, 3)]
  if(fineTune[minIndex] < givenLoss){
    finalModel <- fitKDSN(y = y, X = X, levels = levels1, Dim = Dim1, 
                          sigma = sigma1, lambda = lambda1, alpha = rep(0, levels1), 
                          info = FALSE, 
                          seedW = seedGuess [minIndex, ], 
                          standX = TRUE)
    
    # Include GCV score as attribute
    attr(finalModel, which="GCV") <- fineTune[minIndex]
    return(finalModel)
  }
  else{
    return(estKDSN)
  }
}

# With cross validation
fineTuneCvKDSN <- function(estKDSN, y, X, fineTuneIt=100, info=TRUE, 
                         seedInit=NULL, ncpus=1,
                         cvIndex, lossFunc=devStandard) {
  # Input adjustments
  if(ncpus<1 | floor(ncpus)!=ncpus) {stop("Please supply a positive number of available cpus greater 
                                          or equal to 1 (integer scalar)")}
  fineTune <- vector("numeric", fineTuneIt)
  levels <- estKDSN$Input$levels
  x_new <- c(matrix(c(estKDSN$Input$Dim, estKDSN$Input$sigma, estKDSN$Input$lambda), byrow=TRUE, nrow=3))
  seedOriginal <- estKDSN$Input$seedW
  if(is.null(seedOriginal)) {
    seedOriginal <- seq(0, (levels-1), 1)
  }
  
  # Reproduceability is ensured with seed generation
  set.seed(seedInit)
  seedGuess <- matrix(sample.int(.Machine$integer.max, size=fineTuneIt * levels) * sample(c(-1, 1), size=fineTuneIt * levels, replace=TRUE), nrow=fineTuneIt, ncol=levels)
  envir1 <- environment()
  
  if(ncpus==1) {
    # Fine tune random fourier transformation weights
    for(i in 1:fineTuneIt) {
      fineTune [i] <- lossCvKDSN (parOpt=x_new, y=y, X=X, seedW=seedGuess [i, ],
                                cvIndex=cvIndex, lossFunc=lossFunc)
      if(info) {cat("tuneKDSN", "FineTune =", i, "\n")}
    }
  }
  else {
    tempFunc <- function(i) {
      fineTuneResult <- lossCvKDSN (parOpt=x_new, y=y, X=X, seedW=seedGuess [i, ],
                                    cvIndex=cvIndex, lossFunc=lossFunc)
      return(fineTuneResult)
    }
    
    # Parallel execution with parallel package
    cl <- makeCluster(ncpus)
    clusterExport(cl=cl, varlist=c("seedGuess", "x_new", "y", "X", "gammaPar"), 
                  envir=envir1)
    clusterEvalQ(cl=cl, library(kernDeepStackNet))
    fineTune <- parSapply(cl=cl, X=1:fineTuneIt, FUN=tempFunc)
    stopCluster(cl)
  }
  minIndex <- which.min(fineTune)
  
  # Evaluate loss of given model
  givenLoss <- c(lossCvKDSN (parOpt=x_new, y=y, X=X, cvIndex=cvIndex, seedW=seedOriginal, lossFunc=lossFunc))

  # Refit best model and output if it is better than the given model
  lenx_new <- length(x_new)
  stopifnot((lenx_new%%3) == 0)
  levels1 <- lenx_new/3
  stopifnot(length(seedGuess[minIndex, ]) == levels1)
  Dim1 <- round(x_new[seq(1, lenx_new, 3)])
  sigma1 <- x_new[seq(2, lenx_new, 3)]
  lambda1 <- x_new[seq(3, lenx_new, 3)]
  if(fineTune[minIndex] < givenLoss){
    finalModel <- fitKDSN(y = y, X = X, levels = levels1, Dim = Dim1, 
                          sigma = sigma1, lambda = lambda1, alpha = rep(0, levels1), 
                          info = FALSE, 
                          seedW = seedGuess [minIndex, ], 
                          standX = TRUE)
    # Include GCV score as attribute
    attr(finalModel, which="GCV") <- fineTune[minIndex]
    return(finalModel)
  }
  else{
    return(estKDSN)
  }
}
